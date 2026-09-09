"""Read published Al trajectories; test an overdamped reduction, never invent a clock.

The open dynasor/Zenodo trajectory is a periodic 300 K FCC Al99 NVT crystal,
NOT the production reduced two-row cell and NOT an isolated rigid interface.
We extract adjacent (111) plane-mean displacement differences and examine
their actual correlation. Negative resolved correlation eigenvalues contradict
a reversible harmonic overdamped Markov model even before fitting mobilities.
This research diagnostic cannot emit a production PhysicalTimeCalibration.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import contextlib

import numpy as np


RECORD = 'https://zenodo.org/records/10014454'
TRAJECTORY_URL = RECORD+'/files/dumpT300.NVT.atom.velocity.gz?download=1'
SOURCE_MD5 = 'edf03c7acfd58b8e72ce81db2e6122c0'
MAIN_MD5 = '25f86114b8179f7fcc590b7545acdfa0'


def lammps_frames(stream, maximum_frames):
    """Strict scaled-coordinate orthogonal periodic dump reader; no unit guesses."""
    if maximum_frames < 1:
        raise ValueError('positive frame count required')
    for _ in range(maximum_frames):
        line = stream.readline()
        if not line:
            return
        if line.strip() != 'ITEM: TIMESTEP':
            raise ValueError('expected LAMMPS TIMESTEP')
        step = int(stream.readline())
        if stream.readline().strip() != 'ITEM: NUMBER OF ATOMS':
            raise ValueError('expected atom count')
        count = int(stream.readline())
        if stream.readline().strip() != 'ITEM: BOX BOUNDS pp pp pp':
            raise ValueError('only explicitly orthogonal periodic cells are supported')
        bounds = np.array([list(map(float, stream.readline().split())) for _ in range(3)])
        if bounds.shape != (3, 2) or np.any(~np.isfinite(bounds)) or np.any(np.diff(bounds) <= 0):
            raise ValueError('invalid orthogonal box')
        header = stream.readline().strip().split()
        if header[:2] != ['ITEM:', 'ATOMS']:
            raise ValueError('expected atom columns')
        columns = header[2:]
        if not all(k in columns for k in ('id', 'xs', 'ys', 'zs')):
            raise ValueError('explicit IDs and scaled positions are required')
        raw = ''.join(stream.readline() for _ in range(count))
        values = np.fromstring(raw, sep=' ')
        if values.size != count*len(columns) or np.any(~np.isfinite(values)):
            raise ValueError('truncated/nonfinite atom frame')
        values = values.reshape(count, len(columns))
        ids = values[:, columns.index('id')]
        if np.any(ids != np.rint(ids)) or len(np.unique(ids)) != count:
            raise ValueError('unique integer atom IDs required')
        order = np.argsort(ids)
        positions = values[order][:, [columns.index(k) for k in ('xs', 'ys', 'zs')]]
        yield step, ids[order].astype(np.int64), positions, bounds


class FCCPlaneProjection:
    """A declared periodic cubic crystal, with fixed atom-to-site membership.

    There are n conventional-cell repeats but n (111) plane classes on the
    cubic torus, each containing 4 n^2 atoms. This is not a finite free surface.
    Each q_l is the difference of two plane mean displacements in (normal,
    direct110,in-plane transverse) directions. We DO NOT treat their 4 n^2
    atoms as independent statistical regions or insert an A_c normalization.
    """
    def __init__(self, ids, scaled_positions, bounds, *, repeats, lattice_angstrom):
        if (int(repeats) != repeats or repeats < 2 or not np.isfinite(lattice_angstrom)
                or lattice_angstrom <= 0):
            raise ValueError('physical positive lattice and integer repeats required')
        self.repeats = int(repeats)
        self.lattice_angstrom = float(lattice_angstrom)
        self.ids = np.asarray(ids).copy()
        self.bounds = np.asarray(bounds).copy()
        if len(ids) != 4*repeats**3 or not np.allclose(
                np.diff(bounds).ravel(), repeats*lattice_angstrom, rtol=1e-8, atol=0):
            raise ValueError('declared FCC cubic cell does not match data')
        doubled = np.rint(np.asarray(scaled_positions)*2*repeats).astype(int)
        if np.any(doubled.sum(axis=1) % 2):
            raise ValueError('atoms cannot be assigned to the declared FCC sites')
        modulo = doubled % (2*repeats)
        if len(np.unique(modulo, axis=0)) != len(ids):
            raise ValueError('multiple atoms assigned to one FCC site')
        self.sites = doubled/(2*repeats)
        self.plane = (doubled.sum(axis=1)//2) % repeats
        self.counts = np.bincount(self.plane, minlength=repeats)
        if np.any(self.counts != 4*repeats**2):
            raise ValueError('incorrect periodic (111) plane populations')
        self.basis = np.array([[1, 1, 1], [1, -1, 0], [1, 1, -2]], float)
        self.basis /= np.linalg.norm(self.basis, axis=1)[:, None]

    def evaluate(self, ids, scaled_positions, bounds):
        if not np.array_equal(ids, self.ids) or not np.array_equal(bounds, self.bounds):
            raise ValueError('IDs and box must remain fixed; no silent NPT remapping')
        displacement = np.asarray(scaled_positions)-self.sites
        displacement -= np.rint(displacement)
        displacement *= np.diff(self.bounds).ravel()
        # This is a domain-of-validity refusal, not clipping/dislocation repair.
        if np.max(np.linalg.norm(displacement, axis=1)) >= .25*self.lattice_angstrom:
            raise ValueError('site assignment ambiguous: large displacement/diffusion requires explicit unwrapping')
        means = np.stack([np.bincount(self.plane, weights=displacement[:, j],
                                      minlength=self.repeats)/self.counts for j in range(3)], axis=1)
        return (np.roll(means, -1, axis=0)-means) @ self.basis.T * 1e-10


def lag_covariances(coordinates, lags):
    """Ordered C_ij(k)=<dq_i(t+k)dq_j(t)>; average origins and plane classes.

    Plane classes are correlated: pooling them is not an independent-sample
    count. Each coordinate is centered in its declared time window. No linear
    trend is removed, since doing so could conceal nonstationarity.
    """
    q = np.asarray(coordinates, float)
    lags = np.asarray(lags)
    if (q.ndim != 3 or np.any(~np.isfinite(q)) or lags.ndim != 1
            or np.any(lags != np.rint(lags)) or np.any(lags < 0)
            or np.any(lags >= len(q)) or len(lags) < 1):
        raise ValueError('finite (time,plane,coordinate) and valid integer lags required')
    dq = q-q.mean(axis=0, keepdims=True)
    result = []
    for lag in lags.astype(int):
        later, earlier = (dq[lag:], dq[:-lag]) if lag else (dq, dq)
        result.append(np.einsum('tpi,tpj->ij', later, earlier)/(len(later)*q.shape[1]))
    return np.array(result)


def reversible_harmonic_covariance(hessian, mobility, thermal_energy, times):
    """Exact ordered covariance for a declared harmonic overdamped model.

    H and M use matching physical OR reduced units; this function invents no
    conversion. C0=kT H^-1 and covariance whitening gives
    R_C=sqrt(H) M sqrt(H), not sqrt(M) H sqrt(M) in the same coordinate basis.
    The two symmetric rate matrices have the same eigenvalues, not generally
    the same entries. Return C(t)=kT H^-1/2 exp(-R_C t) H^-1/2.
    """
    H=np.asarray(hessian,float);M=np.asarray(mobility,float);t=np.atleast_1d(np.asarray(times,float))
    if (H.ndim!=2 or H.shape[0]!=H.shape[1] or M.shape!=H.shape
            or any(np.any(~np.isfinite(x)) for x in (H,M,t))
            or not np.allclose(H,H.T) or not np.allclose(M,M.T)
            or not np.isfinite(thermal_energy) or thermal_energy<=0 or np.any(t<0)):
        raise ValueError('finite symmetric positive H,M,thermal energy and nonnegative times required')
    h,V=np.linalg.eigh(H)
    if np.min(h)<=0 or np.min(np.linalg.eigvalsh(M))<=0:
        raise ValueError('stable positive Hessian and mobility required')
    root=(V*np.sqrt(h))@V.T;inverse=(V/np.sqrt(h))@V.T
    rates,U=np.linalg.eigh(root@M@root)
    exp_rate=np.einsum('ij,tj,kj->tik',U,np.exp(-t[:,None]*rates),U)
    return thermal_energy*inverse@exp_rate@inverse


def correlation_audit(coordinates, *, frame_seconds, max_lag, blocks):
    """Empirical window/rate diagnostics, NOT a mobility calibration certificate.

    A reversible overdamped harmonic model has whitened C(t)=exp(-R t),
    R symmetric positive, hence every eigenvalue is positive for finite t.
    We retain antisymmetric sampling error and block ranges independently.
    A maximum block difference is an empirical sensitivity, not a confidence
    interval (blocks can still be correlated).
    """
    q = np.asarray(coordinates, float)
    if not np.isfinite(frame_seconds) or frame_seconds <= 0 or blocks < 2:
        raise ValueError('positive seconds spacing and at least two blocks required')
    width = len(q)//blocks
    if max_lag < 1 or width < 4*(max_lag+1):
        raise ValueError('each time block must contain at least four lag windows')
    lags = np.arange(max_lag+1)
    covariance = lag_covariances(q, lags)
    eig, vectors = np.linalg.eigh(covariance[0])
    if np.any(eig <= 0):
        raise ValueError('positive coordinate covariance required')
    white = (vectors/np.sqrt(eig)) @ vectors.T
    block_cov = np.array([lag_covariances(q[i*width:(i+1)*width], lags)
                          for i in range(blocks)])
    normalized = white @ covariance @ white
    normalized_blocks = white @ block_cov @ white
    symmetric = (normalized+normalized.swapaxes(-1, -2))/2
    symmetric_blocks = (normalized_blocks+normalized_blocks.swapaxes(-1, -2))/2
    eigenvalues = np.linalg.eigvalsh(symmetric)
    block_error = np.max(np.linalg.norm(symmetric_blocks-symmetric, ord=2, axis=(-2, -1)), axis=0)
    antisymmetry = np.linalg.norm(normalized-normalized.swapaxes(-1, -2), ord=2, axis=(-2, -1))/2
    floor = block_error+antisymmetry
    negative = eigenvalues[:, 0] < -floor
    mean_shifts = [(white @ (q[i*width:(i+1)*width].mean(axis=(0, 1))-q.mean(axis=(0, 1)))).tolist()
                   for i in range(blocks)]
    # The spatial sum of periodic adjacent-plane differences vanishes
    # identically. Its mean cannot certify stationarity of individual planes.
    # Retain the old diagnostic but additionally resolve each plane's drift.
    plane_shifts = np.array([(q[i*width:(i+1)*width].mean(axis=0)-q.mean(axis=0)) @ white
                             for i in range(blocks)])
    return dict(lags_seconds=lags*frame_seconds, covariance_m2=covariance,
        whitened_eigenvalues=eigenvalues, block_sensitivity=block_error,
        antisymmetric_residual=antisymmetry, empirical_error_floor=floor,
        negative_beyond_empirical_floor=negative, block_mean_shifts_in_rms_units=mean_shifts,
        spatial_mean_shift_is_not_stationarity_evidence=True,
        block_plane_mean_shifts_in_rms_units=plane_shifts,
        maximum_plane_mean_shift_in_rms_units=float(np.max(np.linalg.norm(plane_shifts,axis=-1))),
        blocks=blocks, block_frames=width, unused_block_tail_frames=len(q)-width*blocks,
        production_calibration_available=False,
        interpretation='negative eigenvalue contradicts reversible harmonic overdamped Markov reduction'
        if np.any(negative) else 'no resolved negative mode in sampled lags; this does NOT establish mobility')


class HashedLimitedReader:
    """Bound network use and hash exactly the compressed prefix actually read."""
    def __init__(self, response, maximum_bytes):
        self.response = response; self.maximum_bytes = maximum_bytes
        self.count = 0; self.digest = hashlib.sha256()

    def read(self, size=-1):
        size = min(size if size >= 0 else 65536, self.maximum_bytes-self.count)
        if size <= 0:
            raise ValueError('declared compressed download budget exhausted')
        value = self.response.read(size)
        self.count += len(value); self.digest.update(value)
        return value


class SequentialRangeReader:
    """Retryable, ordered public byte ranges; bounded prefetch, no whole-file download.

    Range identities are checked so a failed transfer cannot splice unrelated
    bytes or silently restart gzip in its middle. SHA256 covers delivered bytes.
    The upstream whole-file checksum remains unverified for a partial read.
    """
    def __init__(self, url, maximum_bytes, *, chunk_bytes=4*1024**2, workers=4):
        self.url=url; self.maximum_bytes=int(maximum_bytes); self.chunk_bytes=chunk_bytes
        self.pool=ThreadPoolExecutor(max_workers=workers); self.pending={}
        self.chunk=0; self.buffer=b''; self.count=0; self.digest=hashlib.sha256()
        self.total_size=None; self.workers=workers
        for i in range(workers):
            self._submit(i)

    def _submit(self, i):
        start=i*self.chunk_bytes
        if start < self.maximum_bytes:
            end=min(start+self.chunk_bytes,self.maximum_bytes)-1
            self.pending[i]=self.pool.submit(self._fetch,start,end)

    def _fetch(self, start, end):
        for attempt in range(3):
            try:
                request=urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{end}'})
                with urllib.request.urlopen(request,timeout=30) as response:
                    interval=response.headers.get('Content-Range','')
                    if response.status!=206 or not interval.startswith(f'bytes {start}-{end}/'):
                        raise ValueError('server did not return the requested byte interval')
                    total=int(interval.split('/')[-1])
                    data=response.read(end-start+2)
                    if len(data)!=end-start+1:
                        raise IOError('incomplete HTTP byte interval')
                    return data,total
            except (TimeoutError,OSError) as exc:
                if attempt==2:
                    raise IOError(f'public byte-range retrieval failed at {start}: {exc}') from exc

    def read(self,size=-1):
        if size<0:
            size=65536
        result=[]; remaining=size
        while remaining:
            if not self.buffer:
                if self.chunk not in self.pending:
                    raise ValueError('declared compressed download budget exhausted')
                self.buffer,total=self.pending.pop(self.chunk).result()
                if self.total_size is not None and self.total_size!=total:
                    raise ValueError('upstream size changed during transfer')
                self.total_size=total
                self._submit(self.chunk+self.workers); self.chunk+=1
                print(f'public gzip prefix ready: {self.chunk*self.chunk_bytes/1024**2:.0f} MiB',flush=True)
            part=self.buffer[:remaining]; self.buffer=self.buffer[len(part):]
            result.append(part);remaining-=len(part)
        value=b''.join(result);self.count+=len(value);self.digest.update(value)
        return value

    def close(self):
        self.pool.shutdown(wait=True,cancel_futures=True)


def download_projected_trajectory(output, *, frames=2048, max_download_mib=512, ranged=False, source_main_path=None):
    """Stream a bounded prefix; never store the 3.5 GB trajectory in the repo.

    main.in explicitly supplies units metal, dt=.005 ps, dump every5 steps,
    FCC4.065 A,12^3 cells,NVT300K with1ps thermostat damping. Its published
    checksum is verified. Source timestep is not our solver's physical clock.
    The full compressed-file MD5 cannot be verified on a prefix; say so.
    """
    output = Path(output)
    if output.exists():
        raise FileExistsError('fresh output directory required')
    if frames < 2 or max_download_mib <= 0:
        raise ValueError('at least two frames and positive download budget required')
    if source_main_path is None:
        with urllib.request.urlopen(RECORD+'/files/main.in?download=1', timeout=30) as response:
            main = response.read()
    else:
        # Reuse an ACTUAL downloaded artifact on transient network failure.
        # The same published checksum still gates every physical unit claim.
        main=Path(source_main_path).read_bytes()
    if hashlib.md5(main).hexdigest() != MAIN_MD5:
        raise ValueError('published simulation input checksum differs; review source before interpreting units')
    output.mkdir(parents=True)
    # Source artifact, not a handwritten alternative simulation input.
    (output/'source_main.in').write_bytes(main)
    times = []; values = []; first = None; projection = None
    connection=(contextlib.closing(SequentialRangeReader(TRAJECTORY_URL,int(max_download_mib*1024**2)))
                if ranged else urllib.request.urlopen(TRAJECTORY_URL, timeout=45))
    with connection as response:
        limited = response if ranged else HashedLimitedReader(response, int(max_download_mib*1024**2))
        with gzip.GzipFile(fileobj=limited) as compressed:
            with io.TextIOWrapper(compressed, encoding='ascii') as stream:
                for step, ids, positions, bounds in lammps_frames(stream, frames):
                    if projection is None:
                        projection = FCCPlaneProjection(ids, positions, bounds,
                                                        repeats=12, lattice_angstrom=4.065)
                        first = step
                    times.append((step-first)*5e-15)
                    values.append(projection.evaluate(ids, positions, bounds))
                    if len(values) % 128 == 0:
                        np.savez_compressed(output/'projection_checkpoint.npz',
                                            time_seconds=times, coordinates_m=values)
                        print(f'published trajectory frames={len(values)}, compressed MiB={limited.count/1024**2:.1f}', flush=True)
    if len(values) != frames or not np.allclose(np.diff(times), 25e-15, rtol=1e-10, atol=0):
        raise ValueError('incomplete/nonuniform declared trajectory sample')
    np.savez_compressed(output/'plane_coordinates.npz', time_seconds=times, coordinates_m=values)
    metadata = dict(doi='10.5281/zenodo.10014454', source_url=TRAJECTORY_URL,
        creators=['Erik Fransson', 'Paul Erhart'], temperature_K=300,
        thermostat='Nose-Hoover NVT; damping1ps', lattice_angstrom=4.065,
        repeats=12, atoms=6912, atoms_per_plane=576, plane_classes=12,
        timestep_seconds=5e-15, frame_seconds=25e-15, frames=frames,
        duration_seconds=times[-1], first_source_timestep=first,
        main_sha256=hashlib.sha256(main).hexdigest(), main_md5_verified=True,
        trajectory_published_full_md5=SOURCE_MD5, trajectory_full_checksum_verified=False,
        compressed_prefix_bytes=limited.count, compressed_prefix_sha256=limited.digest.hexdigest(),
        http_range_transfer=ranged,
        coordinates='adjacent (111) plane-mean displacement differences: normal,direct110,transverse',
        boundary_conditions='periodic cubic crystal; all atoms mobile; not rigid half-crystal',
        coordinate_equivalent_to_production=False, production_physical_clock=None)
    (output/'source_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n', encoding='utf-8')
    return np.asarray(times), np.asarray(values), metadata

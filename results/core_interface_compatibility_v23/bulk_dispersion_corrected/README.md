The numeric matrices are unchanged. The original `used_in_fit=False` column
was too strong because overlap with earlier selected bulk targets was not
audited. It is corrected to `held_out_certification=False`. This comparison
does not claim blind validation. `bulk_dispersion_audited_partition` is the
complete actual rerun with this explicit provenance qualification.

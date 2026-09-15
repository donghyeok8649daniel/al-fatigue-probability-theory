# AFT 1 surface-load language

Pipeline: selected mesh faces -> stored loads -> AFT text -> parsed IR ->
unit/mesh/tensor validation -> user-approved assignment. Analyzer/chat consumes
a reviewed summary, not arbitrary files or executable commands. This is a
declarative surface-load language, not a completed FVM/FEM solver language.

Example (replace HASH and face indices with those of the actual mesh):

```text
AFT 1
mesh HASH
units mm MPa model_time
load upper
faces 10,11,12
frequency 25.0
parameters 100.0,20.0,0.0,0.0
stress 0,0,0;0,0,0;0,0,normal_mean+normal_amp*sin(2*pi*f*t)
end
balance 1,2,3
```

`parameters` order: normal mean, normal amplitude, shear mean, shear amplitude
(all MPa). Rows/columns are global XYZ. `faces` are zero-based triangle IDs.
Full-line # comments are allowed. Load names must be unique identifiers.
`balance` supplies proposed correction faces; it does not approve correction.
All lengths use mm; f is cycles/model time. Unknown commands, duplicate fields,
missing end, different meshes, asymmetric expressions and unsupported units
are rejected. Imported operator approval and PDE execution are never implicit.

Existing JSON remains supported and uses the same `surface_setup.decode`
validator. AFT is emitted from the same load snapshots, not a second physics
configuration. Solver settings, physical kinetics and mesh generation commands
are deliberately outside v1; no unimplemented `solve` command is accepted.

The local analyzer is deterministic and offline. Its `validated` status means
syntactic surface setup only, never material/kinetic/convergence certification.
External AI answers are advisory drafts. They must go through the same parser
and explicit user review before changing setup; no AI code/tool execution.

## Desktop connection

Face Load -> AFT setup language / Analyzer -> AI Analyzer / Chat.
The editor can read current loads, validate, and explicitly apply reviewed
loads. It does not execute the PDE. The chat window supports separate
conversations and bounded daemon request workers; Tk renders through poll on
the UI thread. Start a new conversation to choose a different explicit API
model ID. Set OPENAI_API_KEY outside the app; keys are not stored in setup.

Setup text, result summary and prior conversation messages are opt-in. The
exact outgoing request is shown before Send. Results are summarized from the
existing local run (including numerical floor when present), never inferred
from the newly edited setup. There are no automatic network requests. Live
API authentication/model access has not been tested; integration uses mocked
transport. store:false is not a claim of zero provider retention. Costs may
apply. Cancelling reception cannot undo an already-sent request.

AI answers remain text. To use a suggested AFT draft, the user must copy it
into the editor, compile it, review it and approve application. Session
histories are in memory only and isolated; they are not Codex conversation
history or shared credentials. No automatic material/kinetic certification.

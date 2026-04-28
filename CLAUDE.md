# Claude Instructions

This project uses **mote** for local issue tracking and agent coordination.

Use mote for:

- work items and follow-up tasks,
- path reservations before editing,
- progress, decision, blocker, and handoff notes,
- completion records.

Do not create new beads (`bd`) issues for project work. Beads are legacy
history only; use `bd` only when explicitly asked to inspect or migrate old
records.

## Basic Workflow

```bash
mote doctor
mote actor show
mote ready
mote show <id>
mote preflight --issue <id> --paths <path> [<path> ...]
mote begin <id> --paths <path> [<path> ...] --note "starting"
mote note <id> --kind progress "what changed"
mote done <id> --note "finished"
```

Use narrow path reservations. If a reservation conflict appears, inspect it with
`mote who-has <path>` and coordinate before editing.

Keep `.mote/` local and out of Git.

## Session Completion

Before stopping, make sure mote items are updated or closed, reservations are
released, intentional tracked changes are committed, and Git is pushed:

```bash
git pull --rebase
git push
git status
```

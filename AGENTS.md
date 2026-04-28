# Agent Instructions

This project uses **mote** for local issue tracking and agent coordination.
Use mote for all new work items, claims, path reservations, notes, handoffs,
and completion records.

## Mote Quick Reference

```bash
mote doctor
mote actor show
mote board
mote ready
mote show <id>
mote preflight --issue <id> --paths <path> [<path> ...]
mote begin <id> --paths <path> [<path> ...] --note "starting"
mote note <id> --kind progress "what changed"
mote done <id> --note "finished"
```

## Mote Rules

- Use `mote` for task tracking and coordination. Do not create new beads
  issues for project work.
- Before editing files, run `mote preflight` for the intended paths and then
  `mote begin` with a narrow path reservation.
- If `mote preflight` or `mote begin` reports a path conflict, inspect it with
  `mote who-has <path>` and coordinate before editing.
- Record material decisions and blockers with `mote note`.
- Finish completed work with `mote done`, or use `mote handoff` / `mote release`
  when stopping before completion.
- Keep `.mote/` local and out of Git. It is intentionally ignored.

## Legacy Beads

This repository previously used **bd** / beads. Active beads have been migrated
to mote; beads are now legacy history only. Do not use `bd` for new tracking
unless explicitly asked to inspect or migrate old bead records.

## Non-Interactive Shell Commands

**ALWAYS use non-interactive flags** with file operations to avoid hanging on
confirmation prompts.

Shell commands like `cp`, `mv`, and `rm` may be aliased to include `-i`
(interactive) mode on some systems, causing the agent to hang indefinitely
waiting for y/n input.

**Use these forms instead:**

```bash
# Force overwrite without prompting
cp -f source dest           # NOT: cp source dest
mv -f source dest           # NOT: mv source dest
rm -f file                  # NOT: rm file

# For recursive operations
rm -rf directory            # NOT: rm -r directory
cp -rf source dest          # NOT: cp -r source dest
```

**Other commands that may prompt:**

- `scp` - use `-o BatchMode=yes` for non-interactive
- `ssh` - use `-o BatchMode=yes` to fail instead of prompting
- `apt-get` - use `-y` flag
- `brew` - use `HOMEBREW_NO_AUTO_UPDATE=1` env var

## Session Completion

Before ending a work session:

1. Record follow-up work in mote.
2. Run relevant quality gates when code changed.
3. Close or update mote items and release reservations.
4. Commit intentional tracked changes.
5. Push Git changes:

```bash
git pull --rebase
git push
git status
```

`git status` must show the branch up to date with `origin` and a clean working
tree before calling the work complete.

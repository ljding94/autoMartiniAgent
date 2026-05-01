# autoMartiniAgent

Agent-agnostic tooling for **automating the all-atom → Martini 3 coarse-grained mapping** problem. Given an atomistic simulation trajectory and a chemical structure, the agent produces a Martini 3 mapping (atom → bead, bead type, bead size) that satisfies Martini sizing rules and yields Gaussian-shaped bond and angle distributions when the AA trajectory is projected through the proposed mapping.

The deliverable is three portable layers:

1. **MCP server** (`mcp_server/`) — exposes the mapping/scoring/repair tools over Model Context Protocol. Mountable in any MCP-aware agent runtime (Claude Code, Codex CLI, OpenCode, Cursor, Continue, …).
2. **Portable skill** (`skill/`) — `SKILL.md` (Claude Code) + `AGENTS.md` (cross-agent community standard); reasoning prose for *when* and *how* to invoke the tools.
3. **Autoresearch program** (`program.md`) — a Karpathy-style agent-runnable protocol that orchestrates the tools end-to-end.

## Status

Early planning. See [`PROGRESS.md`](PROGRESS.md) for the full plan, scope, phases, and status log. See [`program.md`](program.md) for the first-cut autoresearch protocol.

## Background

Triggered by an April 2026 ORNL email chain in which Seonghan Kim circulated a reproducible SMILES → Martini 3 pipeline using [Auto-MartiniM3](https://github.com/Martini-Force-Field-Initiative/Automartini_M3), and Chris Walker stress-tested it on charged polymer monomers (PMETAC, PSBMA) — surfacing two failure modes (subprocess stalls on zwitterionic chemistry, and Martini-3-rule-violating bead sizes) that motivate an agent-driven QA + repair loop on top.

## License

MIT.

# An AI Assistant for Combine

This tutorial introduces an **AI assistant for Combine**: a chat-based agent
that answers Combine questions *with citations* to the real documentation,
paper, source code, and forum, and that can *run* Combine commands to
reproduce, diagnose, and confirm results.

It is meant to lower the barrier to using Combine. You can ask, in plain
language, "how do I set an upper limit on my signal strength?" or "why am I
getting this error from `FitDiagnostics`?", while keeping the answers grounded
in Combine's own sources rather than in a model's memory.

!!! warning "Experimental"
    This is a research prototype, not (yet) an official CMS tool. Always verify what
    it tells you against the linked sources, and review any command before you
    trust its output, especially for anything that affects a physics result.

## How it is built

You do not need to understand the internals to use it, but a quick picture
helps.

<!-- TODO(image): architecture diagram — agent in the middle, retrieval MCP and
     execution MCP on one side, the skill on the other. Add as architecture.png -->
![Architecture of the Combine AI assistant](architecture.png)

These are the main ingredients:

- **The agent**: the chat client you talk to. Two are supported:
  **combagent** (a lightweight terminal client, published on CVMFS so there is
  nothing to install on lxplus) and **Claude Code** (Anthropic's client, handy
  if you already use it locally).
- **MCP servers**: the "tools" the agent is allowed to call. There are two: a
  **retrieval** server that searches five Combine sources (the documentation,
  the [Combine paper](https://arxiv.org/abs/2404.06614), the source code, the
  [cms-talk statistics forum](https://cms-talk.web.cern.ch/c/physics/cat/cat-stats/279),
  and the archived HyperNews forums that preceded it), and an **execution**
  server that runs Combine commands.
- **A skill**: a short set of instructions that tells the agent *how* to use
  Combine: which source to consult for which kind of question, how to read the
  output, and to always cite what it found.

The retrieval server means the agent does not guess: when you ask a question, it
searches the real sources and answers with links you can click to verify.

## Running Combine: two paths

When you ask the agent to actually *run* something, it uses one of two paths, in
order of preference.

1. **Your own Combine environment (recommended).** If you have `combine` on your
   `PATH` (i.e. you sourced a Combine environment *before* launching the agent) it simply runs Combine in your
   current working directory, against your real datacards, with no size limits.
   This is the main mode and the one this tutorial uses.

2. **A remote execution server (fallback).** If Combine is *not* on your `PATH`,
   the agent falls back to a shared execution server running on CERN's cloud
   ([PaaS](https://paas.docs.cern.ch/)). It ships Combine itself and runs your command in an isolated sandbox,
   so it works even with no local Combine install. However, inputs are size-capped
   and timeouts are shorter, so it is best for quick checks rather than large or
   long fits.

!!! tip
    For any real analysis work, use path 1: source your Combine environment
    first. You then get your real files, your real outputs, and no upload
    limits.

## Choosing a model

The agent needs a large language model behind it. Several are pre-configured;
you pick one with the `--model <provider>/<model>` flag (or use the default).

| Provider | How to access | Notes |
|---|---|---|
| **`aigw`** (default) | Your own key from the **CERN AI Gateway** | Recommended for CERN users: data stays within CERN's governed gateway. Serves `qwen3.8-27b-fp16` (the default — CERN-hosted, free to use) and OpenAI's `gpt-5.6-sol/terra/luna-preview` (larger and stronger, but metered). |
| `anthropic` | Your own Anthropic API key | Claude models via the public API. |
| `nrp` | A free token from the [National Research Platform](https://nrp.ai/llms/) | Open-weight models, free for researchers from American universities. |
| `cern-vm` | Nothing (self-hosted) | No key at all, but **CPU-only and very slow** — a last-resort fallback. |

The `litellm` provider, which used to be the default, pointed at the CERN
LiteLLM gateway; that gateway has been **decommissioned** and its models no
longer answer.

For most CERN users the **`aigw`** provider is the right choice. Keys are
per-user (there is no shared key) and must belong to the `cms-combine-agent`
team, so getting one takes four steps:

1. **Subscribe to the e-group**
   [`cms-combine-agent-users`](https://groups-portal.web.cern.ch/group/cms-combine-agent-users/details).

    ![Subscribing to the cms-combine-agent-users e-group](egroup_subscribe.png)

2. **Log in once** at [aigw.cern.ch](https://aigw.cern.ch). The gateway
   registers you on that first visit and adds you to the team a few minutes
   later. Until then you are in the default *sandbox* team, which offers fewer
   models than `cms-combine-agent`.
3. **Create a key** at
   [aigw.cern.ch/ui/api-keys](https://aigw.cern.ch/ui/api-keys/), selecting the
   **`cms-combine-agent`** team. If the team is not offered, step 2 has not gone
   through yet — wait a few minutes and reload. Leave the model selection at
   **All Team Models**.

    ![Creating an API key for the cms-combine-agent team](aigw_api_key.png)

4. **Export it** in the shell where you launch the agent:

    ```shell
    export AIGW_API_KEY=<your-key>
    ```

Putting that `export` in your shell profile saves repeating it every session.
The setup script in step 2 below does **not** fetch a key for you: if
`AIGW_API_KEY` is unset it prints these instructions and nothing else.

To switch model for a session, pass `--model`:

```shell
combagent --model aigw/gpt-5.6-luna-preview
```

## Example 1: combagent on lxplus (main workflow)

This is the recommended way to use the assistant. Everything is on CVMFS, so
there is nothing to install.

**Step 1: source a Combine environment** so that `combine` is on your `PATH`.
Use your own CMSSW area with Combine built in it (see the
[installation instructions](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#within-cmssw-recommended-for-cms-users)):

```shell
cd /path/to/CMSSW_16_0_0/src
cmsenv
cd /path/to/your/analysis    # where your datacards live
```

**Step 2: source the assistant setup.** This puts the `combagent` client (our
maintained fork of opencode) on your `PATH` and wires in the Combine tools,
skill, and model configuration:

```shell
source /cvmfs/cms-griddata.cern.ch/cat/sw/combine-assistant/latest/bin/setup.sh
export AIGW_API_KEY=<your-key>
```

The key is yours to create — see [Choosing a model](#choosing-a-model) above —
and the setup script does not fetch one for you. If `AIGW_API_KEY` is unset it
prints those same instructions instead.

**Step 3: launch the agent:**

```shell
combagent
```

You should see the **combagent** banner, and you can start chatting.

<!-- TODO(image): screenshot of the combagent banner / a first prompt.
     Add as combagent_lxplus.png -->
![combagent running on lxplus](combagent_banner.png)

Some things to try:

- *"What does the `--robustFit` option do, and when should I use it?"* — the
  agent searches the docs and answers with a link.
- *"Run an asymptotic upper limit on `datacard.txt` and report the expected
  limit."* — because you sourced Combine in step 1, it runs
  `combine -M AsymptoticLimits datacard.txt` in your directory and reports
  the result.
- *"I get `Error: ... covariance matrix ...` from FitDiagnostics — what's going
  on?"* — the agent checks the forum and code and proposes a fix.

!!! note
    The `aigw` models and the self-hosted `cern-vm` require the CERN network
    (lxplus/SWAN, or VPN).

## Example 2: Claude Code, locally

If you already use [Claude Code](https://claude.com/claude-code) on your own
machine, you can use the same Combine tools and skill from there.

**Step 1: get the assistant configuration** (it ships the tool registrations
and the skill):

```shell
git clone https://github.com/maxgalli/combine-assistant.git
cd combine-assistant
```

**Step 2: (optional) source a Combine environment** in the same shell if you
want commands to run against a local Combine install; otherwise the remote
execution server is used automatically.

**Step 3: launch Claude Code** in that directory:

```shell
claude
```

Claude Code automatically registers the Combine retrieval and execution servers
and discovers the skill, and it uses its own (Anthropic) model, so there is no
separate model key to configure. You can then ask the same kinds of questions as
above.

!!! tip
    On a laptop with no Combine installed and no CERN network, the agent can
    still *answer* questions (retrieval works over the public sources) but
    cannot *run* Combine. Source a Combine environment, or use it from lxplus,
    to enable execution.

## Where to go next

- Ask the agent itself how to do a specific task (that is what it is for).
- Cross-check its citations: every answer should point at the docs, paper,
  code, or forum. If it cannot find something, it will say so rather than
  guess.
- For the underlying Combine methods it drives, see the rest of these
  tutorials and the [documentation home page](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/).

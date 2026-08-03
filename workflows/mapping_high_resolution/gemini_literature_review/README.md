# Gemini literature-review client

This dedicated folder supports source-grounded literature work for Manuscript 3
without storing credentials in the repository.

## Authentication

1. Create an API key in Google AI Studio: https://aistudio.google.com/apikey
2. In Windows, open **Environment Variables**.
3. Under **User variables**, create `GEMINI_API_KEY` and paste the key as its
   value.
4. Close and reopen VS Code/Codex so the new environment variable is loaded.

Do not paste the key into a script, chat, committed file, or `.env.example`.
The local `.env` filename is ignored as an additional safeguard, but the
scripts intentionally use the Windows environment variable.

## Connection test

From this folder, run:

```powershell
npm run check
```

The default model can be changed through a `GEMINI_MODEL` environment variable.
Generated literature-review outputs belong in `outputs/`; this directory is
ignored except for its placeholder.

## Design decisions

- Official `@google/genai` SDK.
- Stable Gemini API (`v1`).
- Interactions API for new development.
- No API key in source code.
- Dedicated scripts and outputs so the analysis repository is not cluttered.

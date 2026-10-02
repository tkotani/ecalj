# ecaljdoc

## What is ecaljdoc?

ecaljdoc is a suite of documents for the [ecalj package](https://github.com/tkotani/ecalj/).
Since we use VitePress for ecaljdoc, we have:

* Read ecaljdoc — <https://ecalj.github.io/ecaljdoc>
* Source ecaljdoc — <https://github.com/ecalj/ecaljdoc>

If you can join to edit ecaljdoc, ask msobt at <https://github.com/msobt>.

> 🤖 **Auto-edited sections (2026-05)** — the TOML-migration banners on
> manual / install / implementation pages and the new
> [`manual/toml_migration.md`](./manual/toml_migration.md) hub were
> drafted by Claude (Opus 4.7) in a single Claude Code session
> (commit `86e7646`). The factual structure (file roles, Legacy2toml.py
> usage, migrated sample tree) reflects the actual ecalj state, but
> commit hashes, line counts, and English wording may need a human
> review pass. Treat with care; correct as you find issues.

---

### How to edit this document with VitePress

* Get the source of ecaljdoc:

  ```bash
  git clone git@github.com:ecalj/ecaljdoc.git
  ```

  Changes at `git@github.com:ecalj/ecaljdoc.git` are promptly reflected at <https://ecalj.github.io/ecaljdoc>.

* Install the dependencies:

  ```bash
  npm install
  ```

* Generate the document site (ecaljdoc) on your local machine:

  ```bash
  npm run docs:dev
  ```

  If it works, you will see the following message:

  ```text
  vitepress v1.6.3

  ➜  Local:   http://localhost:5173/ecaljdoc/
  ➜  Network: use --host to expose
  ➜  press h to show help
  ```

  The generated document site can be seen at `http://localhost:5173/ecaljdoc/`. Open the URL in your web browser.

Note that we have two key files for VitePress:

* Main page configuration — `index.md`
* Sidebar configuration — `.vitepress/config.mts`

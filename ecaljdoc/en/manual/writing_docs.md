# Writing ecaljdoc {#writing-docs}

> English version of the page [manual/writing_docs](../../manual/writing_docs.md). The pages at the top of ecaljdoc are the reference; this one was made from manual/writing_docs.md of commit e1dd775c3 (source blob 79c5fd0d8e8e).

Notes for writing and editing this site (ecaljdoc). They hold for people and for AI alike.

## Where it lives, and publishing {#where-and-publishing}

- The reference is `ecaljdoc/` in the ecalj repository. Code and documents may be changed in the same commit (commit messages in English).
- The public site (https://ecalj.github.io/ecaljdoc/) is built by GitHub Pages from the GitHub repository `ecalj/ecaljdoc`. To publish, the committed
  `ecaljdoc/` of the ecalj HEAD is copied there with `TOOLS/publish_ecaljdoc.sh` (`--push` to push). Edits made directly in `ecalj/ecaljdoc` are overwritten by the next publication.
- Three things are published in order: the ecalj push, ecaljdoc, and the QSGW band database; `TOOLS/publish_all.sh` does all. Push only when the maintainer says so.

## Preview locally {#preview}

```bash
cd ecalj/ecaljdoc
npm install            # first time only
npm run docs:dev       # http://localhost:5173/ecaljdoc/
npx vitepress build    # check that the site builds
```

## Pages and the menu {#pages-and-menu}

- Every page is in the menu on the left. The menu is the only table of contents; no page keeps a second copy of that list.
- Each piece of content is written in one place. Where the text touches the subject of another page, link to that page (to its section if possible); do not copy the content.
- When you add a page, put it into `sidebarRoot` in `.vitepress/config.mts`. If the item name is Japanese, add its English name to `enText`.

## Japanese and English {#japanese-and-english}

- The Japanese pages are the reference. The English version is the page of the same name under `en/` (switch with the button at the top right).
- When you change a Japanese page, change the English one too. Line 3 of an English page carries the mark of its source (source blob), and
  `python3 TOOLS/ecaljdoc_en.py status` lists the English pages that are behind. When it prints nothing, they agree. After editing, update that mark.
- For a new page, `python3 TOOLS/ecaljdoc_en.py copy <page>.md` makes the English skeleton; then translate the text.

## Style {#style}

- Half of the readers are AI, but write so that people can read it: balanced and short, without long disclaimers or repetitions meant for AI.
- Write only the reasons and the cautions that are needed. How we got there (the history) does not go here (it goes to the development records).
- Equations carry numbers (`\tag{n}`) and are referred to by number; figures and tables are numbered too.
- Do not separate who wrote what (a person or an AI) or who asked for it; git tells that.
- In the Japanese pages, core is written コア (not 芯); semicore is セミコア and core hole is コアホール, fractional Z is フラクショナルな Z (not 半端な Z). Some other terms are also written in katakana.
- Section numbers ("1." for h2, "1.2" for h3) are added by CSS (`.vitepress/theme/custom.css`), and the outline on the right shows the same numbers; use them in discussions. Do not write numbers in the headings by hand. Refer to a section in the text by its automatic number ("§3.2"); when sections are added or reordered, fix the "§" references to that page (`pageClass: secnum-manual` in the frontmatter stops the automatic numbers, but no page uses it now). Avoid headings that begin with a numeral (they read as part of the number): write "three-step procedure" rather than "3-step procedure". Use only one h1, the title of the page (a second one shifts the numbers). Give a heading made only of math a `{#id}` (its name would be empty, it would drop out of the outline and shift the numbers).
- For an input key whose template value (written by `ctrlgenToml.py`) differs from the default of the code (used when the key is absent), write both values.

## VitePress pitfalls {#vitepress-pitfalls}

- `<...>` outside code marks (backquotes) in the text stops the build; write `&lt;` or put it in code.
- A double curly brace (two `{` in a row) stops the build even inside inline code marks (Vue takes it for an expression); inside a code block (```) it is fine.
- `<` inside inline math stops the build; use `\lt`. A space before the closing `$` stops it too.
- A bare URL inside Japanese text takes the following characters into the link; write `[text](URL)`.
- In menu item names write `<sname>` as `&lt;sname&gt;`.

## Links {#links}

- From a page to outside ecaljdoc (Samples, sources), use GitHub URLs (relative `../` links break on the site). Do not link to the development records (`MD/` of the repository).
- Check with `python3 TOOLS/doclinks.py`.

## Visible corrections (red pen) {#red-pen}

When corrections are to be checked, make them visible: the old text in `<del class="clfix">…</del>`, the new text in `<span class="clfix">…</span>`.
For a whole equation or paragraph, wrap the old one in `<div class="clfix-old">` and the new one in `<div class="clfix">` (blank lines inside). After the check, remove the marks.

## Removed things {#removed}

Pages and code are not deleted: they are moved to the local `trash/` and untracked. For code, add a line with the commit before removal to [Removed code](./past_code).

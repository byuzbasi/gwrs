# Build the documentation locally

The website covers joint subset/kernel/bandwidth selection for ordinary GWR,
GWR-Ridge, GWR-Lasso, GWR-ENET, GWR-SCAD and GWR-MCP, their supporting
spatial inputs, diagnostics, prediction and the documented county data.
`build-site.R` selects the website reference topics without modifying any
package function, help file, default or numerical implementation.

With the current package and the declared documentation dependencies already
installed, run this command from the repository root:

```sh
Rscript --vanilla pkgdown/build-site.R
```

Open `docs/index.html` to inspect the generated site. The build executes the
small five-method and ordinary-GWR selection guides. These are demonstrations,
separate from the paper's simulation or empirical study.
The other articles display commands and existing results without fitting.

Generated pages are local previews. This command does not deploy a website or
start GitHub Actions. Deployment and repository visibility are separate steps.
The site uses local image and theme assets; no hosted badge is needed.

The banner and logo are original vector artwork. Their contours illustrate
the package's spatial theme and do not represent estimated coefficient values.
The county selection image is copied unchanged from the checked manuscript
figure. Its source filename and SHA-256 are recorded in
`figure-provenance.json`; the underlying study is linked from the article.

# Change Log

## [Unreleased](https://github.com/sanger-pathogens/companion/tree/HEAD)

[Full Changelog](https://github.com/sanger-pathogens/companion/compare/v1.0.1...HEAD)

**Fixed bugs:**

- Run fails if Aragorn predicts unreliable anticodons [\#61](https://github.com/sanger-pathogens/companion/issues/61)

**Closed issues:**

- sudo [\#57](https://github.com/sanger-pathogens/companion/issues/57)

**Merged pull requests:**

- Final \(?\) documentation rewrite [\#89](https://github.com/sanger-pathogens/companion/pull/89) ([trstickland](https://github.com/trstickland))
- Doc rewrite [\#88](https://github.com/sanger-pathogens/companion/pull/88) ([trstickland](https://github.com/trstickland))
- travis CI now pulls rather than builds the docker image [\#86](https://github.com/sanger-pathogens/companion/pull/86) ([trstickland](https://github.com/trstickland))
- travis build fixed [\#83](https://github.com/sanger-pathogens/companion/pull/83) ([trstickland](https://github.com/trstickland))
- snap wasn't in PATH; symlink from /usr/local/bin [\#82](https://github.com/sanger-pathogens/companion/pull/82) ([trstickland](https://github.com/trstickland))
- Docker build fail [\#81](https://github.com/sanger-pathogens/companion/pull/81) ([trstickland](https://github.com/trstickland))
- update readme format [\#79](https://github.com/sanger-pathogens/companion/pull/79) ([ssjunnebo](https://github.com/ssjunnebo))
- add AUTHORS file [\#78](https://github.com/sanger-pathogens/companion/pull/78) ([ssjunnebo](https://github.com/ssjunnebo))
- add systematic ID [\#75](https://github.com/sanger-pathogens/companion/pull/75) ([satta](https://github.com/satta))
- add multiple transcript qualifiers [\#73](https://github.com/sanger-pathogens/companion/pull/73) ([satta](https://github.com/satta))
- ensure that reports step always finishes with required files [\#72](https://github.com/sanger-pathogens/companion/pull/72) ([satta](https://github.com/satta))
- ignore nonzero exit codes from speck [\#71](https://github.com/sanger-pathogens/companion/pull/71) ([satta](https://github.com/satta))
- Repo cleanup [\#70](https://github.com/sanger-pathogens/companion/pull/70) ([satta](https://github.com/satta))
- Stability and robustness fixes [\#69](https://github.com/sanger-pathogens/companion/pull/69) ([satta](https://github.com/satta))
- Robustness improvements [\#67](https://github.com/sanger-pathogens/companion/pull/67) ([satta](https://github.com/satta))
- major improvements to EMBL converter [\#66](https://github.com/sanger-pathogens/companion/pull/66) ([satta](https://github.com/satta))
- separate genes and pseudogenes in stats generation [\#65](https://github.com/sanger-pathogens/companion/pull/65) ([satta](https://github.com/satta))
- use correct package name in testing [\#64](https://github.com/sanger-pathogens/companion/pull/64) ([satta](https://github.com/satta))
- EMBL export improvements [\#63](https://github.com/sanger-pathogens/companion/pull/63) ([satta](https://github.com/satta))
- Bug fixes, cleanups and robustness improvements [\#62](https://github.com/sanger-pathogens/companion/pull/62) ([satta](https://github.com/satta))
- Include UTRs in EMBL conversion [\#60](https://github.com/sanger-pathogens/companion/pull/60) ([satta](https://github.com/satta))
- Channel `hints` channel is never used. [\#59](https://github.com/sanger-pathogens/companion/pull/59) ([pditommaso](https://github.com/pditommaso))
- always truncate inputs before sanitizing [\#58](https://github.com/sanger-pathogens/companion/pull/58) ([satta](https://github.com/satta))

## [v1.0.1](https://github.com/sanger-pathogens/companion/tree/v1.0.1) (2016-04-05)
[Full Changelog](https://github.com/sanger-pathogens/companion/compare/v1.0.0...v1.0.1)

**Implemented enhancements:**

- Use preferred product name in EMBL output [\#51](https://github.com/sanger-pathogens/companion/issues/51)
- Download option for table content [\#47](https://github.com/sanger-pathogens/companion/issues/47)

**Fixed bugs:**

- Proactive sanitization of input headers with special characters [\#49](https://github.com/sanger-pathogens/companion/issues/49)

**Merged pull requests:**

- Address problems with EMBL output  [\#56](https://github.com/sanger-pathogens/companion/pull/56) ([satta](https://github.com/satta))
- small fixes [\#54](https://github.com/sanger-pathogens/companion/pull/54) ([satta](https://github.com/satta))
- various improvements [\#53](https://github.com/sanger-pathogens/companion/pull/53) ([satta](https://github.com/satta))
- use new version of libgsl package [\#50](https://github.com/sanger-pathogens/companion/pull/50) ([satta](https://github.com/satta))

## [v1.0.0](https://github.com/sanger-pathogens/companion/tree/v1.0.0) (2016-01-28)
**Implemented enhancements:**

- Allow optional alphanumeric random locus tags [\#45](https://github.com/sanger-pathogens/companion/issues/45)

**Fixed bugs:**

- Missing SVG Perl module [\#1](https://github.com/sanger-pathogens/companion/issues/1)

**Merged pull requests:**

- ENA validation and ID assignment [\#46](https://github.com/sanger-pathogens/companion/pull/46) ([satta](https://github.com/satta))
- 'Finish line' fixes towards ENA validity [\#43](https://github.com/sanger-pathogens/companion/pull/43) ([satta](https://github.com/satta))
- allow pipeline to complete when no genes are annotated at all [\#42](https://github.com/sanger-pathogens/companion/pull/42) ([satta](https://github.com/satta))
- Make sure Docker Hub builds working images [\#41](https://github.com/sanger-pathogens/companion/pull/41) ([satta](https://github.com/satta))
- Stability improvements [\#40](https://github.com/sanger-pathogens/companion/pull/40) ([satta](https://github.com/satta))
- fix Circos drawing [\#39](https://github.com/sanger-pathogens/companion/pull/39) ([satta](https://github.com/satta))
- use whole genome as RATT input, not just chromosomes [\#38](https://github.com/sanger-pathogens/companion/pull/38) ([satta](https://github.com/satta))
- use new Docker hub container [\#37](https://github.com/sanger-pathogens/companion/pull/37) ([satta](https://github.com/satta))
- pseudogene and chromosome handling [\#36](https://github.com/sanger-pathogens/companion/pull/36) ([satta](https://github.com/satta))
- skip Pfam hits with invalid converted ranges [\#35](https://github.com/sanger-pathogens/companion/pull/35) ([satta](https://github.com/satta))
- Latest work [\#34](https://github.com/sanger-pathogens/companion/pull/34) ([satta](https://github.com/satta))
- small fixes [\#33](https://github.com/sanger-pathogens/companion/pull/33) ([satta](https://github.com/satta))
- Robustness improvement [\#32](https://github.com/sanger-pathogens/companion/pull/32) ([satta](https://github.com/satta))
- Latest work [\#31](https://github.com/sanger-pathogens/companion/pull/31) ([satta](https://github.com/satta))
- do not treat lowercase input sequences as repeat masked in LAST [\#30](https://github.com/sanger-pathogens/companion/pull/30) ([satta](https://github.com/satta))
- AGP output and RATT integration [\#29](https://github.com/sanger-pathogens/companion/pull/29) ([satta](https://github.com/satta))
- fixes and improvements [\#28](https://github.com/sanger-pathogens/companion/pull/28) ([satta](https://github.com/satta))
- fixes and improvements [\#27](https://github.com/sanger-pathogens/companion/pull/27) ([satta](https://github.com/satta))
- fixes and improvements [\#26](https://github.com/sanger-pathogens/companion/pull/26) ([satta](https://github.com/satta))
- fixes and improvements [\#25](https://github.com/sanger-pathogens/companion/pull/25) ([satta](https://github.com/satta))
- some polishing of pipeline output and postprocessing [\#24](https://github.com/sanger-pathogens/companion/pull/24) ([satta](https://github.com/satta))
- make AGP splitting seqid comparisons more stringent [\#23](https://github.com/sanger-pathogens/companion/pull/23) ([satta](https://github.com/satta))
- only keep metadata in references.json [\#22](https://github.com/sanger-pathogens/companion/pull/22) ([satta](https://github.com/satta))
- integrate ab initio models into reference directory [\#21](https://github.com/sanger-pathogens/companion/pull/21) ([satta](https://github.com/satta))
- update reference examples [\#20](https://github.com/sanger-pathogens/companion/pull/20) ([satta](https://github.com/satta))
- make sure to fix ownership on Docker [\#19](https://github.com/sanger-pathogens/companion/pull/19) ([satta](https://github.com/satta))
- do not fail on rm [\#18](https://github.com/sanger-pathogens/companion/pull/18) ([satta](https://github.com/satta))
- clean up after RATT [\#17](https://github.com/sanger-pathogens/companion/pull/17) ([satta](https://github.com/satta))
- add features and cleanup [\#16](https://github.com/sanger-pathogens/companion/pull/16) ([satta](https://github.com/satta))
- allow transcript evidence [\#15](https://github.com/sanger-pathogens/companion/pull/15) ([satta](https://github.com/satta))
- small fixes/improvements [\#14](https://github.com/sanger-pathogens/companion/pull/14) ([satta](https://github.com/satta))
- make weight definition for integration step user configurable [\#13](https://github.com/sanger-pathogens/companion/pull/13) ([satta](https://github.com/satta))
- some improvements [\#12](https://github.com/sanger-pathogens/companion/pull/12) ([satta](https://github.com/satta))
- some work on stability and variety [\#11](https://github.com/sanger-pathogens/companion/pull/11) ([satta](https://github.com/satta))
- simplify and improve EMBL export for RATT [\#10](https://github.com/sanger-pathogens/companion/pull/10) ([satta](https://github.com/satta))
- fix Circos bin drawing  [\#9](https://github.com/sanger-pathogens/companion/pull/9) ([satta](https://github.com/satta))
- do not treat reference chromosome 'numbers' as numbers [\#8](https://github.com/sanger-pathogens/companion/pull/8) ([satta](https://github.com/satta))
- add tantan, make PTU smoothing optional [\#7](https://github.com/sanger-pathogens/companion/pull/7) ([satta](https://github.com/satta))
- use correct Travis build badge [\#6](https://github.com/sanger-pathogens/companion/pull/6) ([satta](https://github.com/satta))
- use $baseDir consistently [\#5](https://github.com/sanger-pathogens/companion/pull/5) ([satta](https://github.com/satta))
- Manual merge of pull request \#2 [\#4](https://github.com/sanger-pathogens/companion/pull/4) ([satta](https://github.com/satta))
- Use into operator with closure syntax [\#3](https://github.com/sanger-pathogens/companion/pull/3) ([pditommaso](https://github.com/pditommaso))



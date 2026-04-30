# GxEScanR 2.0.0

* Added a `NEWS.md` file to track changes to the package.

# GxEScanR 2.0.1

* Updated documentation by adding meanings of acronyms.
* Added tests to completely check outputs.

# GxEScanR 2.0.2

* Fixed issue with compiling on solaris.

# GxEScanR 2.0.0.9000

* Development version. Redesigned to depend on the lsReg and BinaryDosage
  packages; C++ source code moved to lsReg.
* Replaced gwas() and gweis() with gweis.mem() and rungweis() to support
  a wider set of GWEIS tests including gene-only, GxE interaction, joint,
  E|G, case-only, and control-only models.

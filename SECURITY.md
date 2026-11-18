# Security policy

## Supported versions

Security fixes go into the latest release and the `main` branch.

## Reporting a vulnerability

Please do not open a public issue. Report it privately through GitHub: on the repository's **Security** tab, choose
**Report a vulnerability** ([GitHub's guide](https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing-information-about-vulnerabilities/privately-reporting-a-security-vulnerability)).
Include the affected version or commit, the input or command that triggers it, and what happens.

You should get an answer within two weeks. Once a fix is released, the advisory is published with credit to the
reporter, unless you ask otherwise.

varcall-rs reads BAM, CRAM, FASTQ and VCF files and runs no network services. Of most interest are inputs that make
it read or write out of bounds, run code, or write somewhere other than its output paths.

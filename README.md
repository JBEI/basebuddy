<p align="center">
  <img src="resources/logo/svg/logo-no-background.svg" />
</p>

A web app for recoding genes by organism.

See the deployment at [basebuddy.lbl.gov](https://basebuddy.lbl.gov).

## Citing This Work

If you use this tool, please cite the BaseBuddy paper (Schmidt et al. 2023) and DnaChisel:

> Zulkower, Valentin, and Susan Rosser. 2019. “DNA Chisel, a Versatile Sequence Optimizer.” bioRxiv. https://doi.org/10.1101/2019.12.16.877480.

## To Run Locally

* Clone the repo: `git clone git@github.com:JBEI/basebuddy.git`
* Decompress the cocoput table: `cd data && tar -xvf cocoput_table.tsv.tar.gz && cd ..`
* Run streamlit in docker: `bash run_locally.sh`

## Acknowledgements

BaseBuddy utilizes many separate libraries and packages including:

- [DNAChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)
- [Streamlit](https://github.com/streamlit/streamlit)

We thank all their contributors and maintainers!

Use of the third-party software, libraries or code BaseBuddy may be governed by separate terms and conditions or license provisions. Your use of the third-party software, libraries or code is subject to any such terms and you should check that you can comply with any applicable restrictions or terms and conditions before use.

## License

BaseBuddy is distributed under a modified BSD license (see LICENSE).

## Copyright Notice

BaseBuddy Copyright (c) 2023, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy) and University 
of California, Berkeley. All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.

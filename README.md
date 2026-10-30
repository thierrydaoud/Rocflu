# Rocflu

The RocfluMP flow solver coupled to the ppiclF particle tracker.

📘 **Documentation: https://thierrydaoud.github.io/Rocflu/**

```bash
./rocflu.sh build                                   # build everything
./rocflu.sh test                                    # regression suite
./rocflu.sh run /path/to/case/<case>.inp -n 4       # check, build, pre-process and submit a case
./rocflu.sh post /path/to/case/<case>.inp           # results -> ParaView files
```

Also: the [RocfluMP manual](docs/Rocflu_manual.pdf) (PDF) and the
[ppiclF documentation](https://dpzwick.github.io/ppiclF-doc/user/external.html).
The previous version of the code is on the [`legacy`](https://github.com/thierrydaoud/Rocflu/tree/legacy) branch.

# KKMCee-dev
KKMCee Development

```bash
# clone the KKMCee-dev git repository and then ...
cd KKMCee-dev
ln -f -s dizet-6.42 dizet # note: use dizet-6.42 or dizet-6.21
nice autoreconf -f -i 2>&1 | tee autoreconf.out.txt
# see "./configure --help" for available options
nice ./configure 2>&1 | tee configure.out.txt
nice make 2>&1 | tee make.out.txt
```

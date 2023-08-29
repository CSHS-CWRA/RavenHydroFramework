# Raven Standalone Examples

The examples assume that you have access to an execuable version of Raven.

## Alouette5_60hours

Model for the two first days of 2016 at *hourly* temporal resolution (hypothetical rainfall).

Single basin segmented into 23 HRUs.

In a bash terminal (Linux):

```bash
$ mkdir -p Alouette5_60hours/output
$ cd Alouette5_60hours/output
$ Raven ../Alouette5 -o . -n
```
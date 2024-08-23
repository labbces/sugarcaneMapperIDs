# sugarcaneMapperIDs

The goal of this app is the following: given a DNA or protein sequence for sugarcane map it to a desired database (from sequences in Conekt Grasses) and return the corresponding identifier used to compute expression profiles or co-expression networks.

## Required software

```bash
mkdir sugarcaneMapperIDs
cd sugarcaneMapperIDs
python3 -m venv .venv
. .venv/bin/activate
pip install flask
pip install peewee
```

## Installing DB schema

```bash
cat db/schema.sql | sqlite3 sugarcaneSequences.db
```

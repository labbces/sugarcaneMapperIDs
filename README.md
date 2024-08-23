# sugarcaneMapperIDs

The goal of this app is to given a DNA or protein sequence for sugarcane map it to a desired database (from sequences in Conekt Grasses) and give the corresponding identifier used to compute expression profiles or co-expression networks.

## Required software
```bash
pip install flask
pip install peewee
```

## Installing DB schema

```bash
cat db/schema.sql | sqlite3 sugarcaneSequences.db
```

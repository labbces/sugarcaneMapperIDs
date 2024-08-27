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

The schema is available as SQL, and can be used to initialized the database in the following way:

```bash
cat db/schema.sql | sqlite3 sugarcaneSequences.db
```

The model.py also knows how to implement the schema, and hte preferred way to initialize the database is the following:

```bash
python3 db/model.py
```

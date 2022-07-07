CREATE TABLE ref (
  id INTEGER PRIMARY KEY,
  name TEXT NOT NULL UNIQUE
);
CREATE TABLE bam (
  id INTEGER PRIMARY KEY,
  uuid TEXT NOT NULL UNIQUE,
  filepath TEXT
);
CREATE TABLE chunk (
  bam_id INTEGER NOT NULL,
  ref_id INTEGER NOT NULL,
  bin_id INTEGER NOT NULL,
  chunk_beg INTEGER NOT NULL,
  chunk_end INTEGER NOT NULL,
  FOREIGN KEY(bam_id) REFERENCES bam(id),
  FOREIGN KEY(ref_id) REFERENCES ref(id)
);
CREATE TABLE intv (
  bam_id INTEGER NOT NULL,
  ref_id INTEGER NOT NULL,
  -- WINDOW is a keyword in postgres and needs quotes as column name
  "window" INTEGER NOT NULL,
  ioffset INTEGER NOT NULL,
  FOREIGN KEY(bam_id) REFERENCES bam(id),
  FOREIGN KEY(ref_id) REFERENCES ref(id)
);
CREATE INDEX chunk_bin_ref_bam_idx ON chunk(bin_id, ref_id, bam_id);
CREATE INDEX intv_window_ref_bam_idx ON intv("window", ref_id, bam_id);

pub const BED3: &str = r#"
table bed3
"Simple bed"
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
)
"#;

pub fn bed_autosql(rest: &str) -> String {
    let extra_fields = if rest.is_empty() {
        0
    } else {
        rest.split('\t').count()
    };
    let mut def = "\
table bed
\"Browser Extensible Data\"
(
    string chrom;       \"Reference sequence chromosome or scaffold\"
    uint   chromStart;  \"Start position in chromosome\"
    uint   chromEnd;    \"End position in chromosome\"
"
    .to_string();
    const FIELDS: &[&str] = &[
        "   string name;        \"Name of item.\"\n",
        "   uint score;          \"Score (0-1000)\"\n",
        "   char[1] strand;     \"+ or - for strand\"\n",
        "   uint thickStart;   \"Start of where display should be thick (start codon)\"\n",
        "   uint thickEnd;     \"End of where display should be thick (stop codon)\"\n",
        "   uint reserved;     \"Used as itemRgb as of 2004-11-22\"\n",
        "   int blockCount;    \"Number of blocks\"\n",
        "   int[blockCount] blockSizes; \"Comma separated list of block sizes\"\n",
        "   int[blockCount] chromStarts; \"Start positions relative to chromStart\"\n",
        "   int expCount;	\"Experiment count\"\n",
        "   int[expCount] expIds;	\"Comma separated list of experiment ids. Always 0,1,2,3....\"\n",
        "   float[expCount] expScores; \"Comma separated list of experiment scores.\"\n",
    ];
    for field in &FIELDS[0..extra_fields.min(FIELDS.len())] {
        def.push_str(field);
    }
    for i in FIELDS.len()..extra_fields.max(FIELDS.len()) {
        def.push_str(&format!(
            "   lstring field{};	\"Undocumented field\"\n",
            i + 3 + 1
        ))
    }
    def.push(')');
    def
}

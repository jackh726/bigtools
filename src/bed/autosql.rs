/*!
Utitilies for reading and writing the autosql section of a bigBed.
*/

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

// Defined by https://github.com/ucscGenomeBrowser/kent/blob/c26640b68ba8ad219e7d79c3f8251ea20f9f57e0/src/hg/autoSql/autoSql.doc
pub mod parse {
    mod parser {
        pub(super) struct Parser<'a> {
            pub(super) data: &'a str,
            pub(super) start_cursor: usize,
            pub(super) end_cursor: usize,
        }

        impl<'a> Parser<'a> {
            pub(super) fn of(data: &'a str) -> Self {
                Parser {
                    data,
                    start_cursor: 0,
                    end_cursor: 0,
                }
            }

            fn is_word_delimiter(data: char) -> bool {
                data.is_whitespace()
                    || data == ';'
                    || data == '('
                    || data == ')'
                    || data == '['
                    || data == ']'
                    || data == ','
            }

            fn take_whitespace(&mut self) {
                loop {
                    let mut chars = self.data[self.start_cursor..].char_indices().peekable();
                    let next_char = match chars.next() {
                        Some(char) => char,
                        None => {
                            self.end_cursor = self.start_cursor;
                            return;
                        }
                    };
                    if !next_char.1.is_whitespace() {
                        break;
                    }
                    match chars.peek() {
                        Some((next_index, _)) => {
                            self.start_cursor = next_index + self.start_cursor;
                            self.end_cursor = self.start_cursor;
                        }
                        None => {
                            self.start_cursor = self.data.len();
                            self.end_cursor = self.start_cursor;
                        }
                    }
                }
            }

            /// Advances the start cursor to the end of the last peeked value
            pub(super) fn take(&mut self) -> &'a str {
                let start = self.start_cursor;
                let end = self.end_cursor;
                self.start_cursor = self.end_cursor;
                &self.data[start..end]
            }

            /// First advances the start cursor past any whitespace.
            /// Then, returns a peek of contiguos non-whitespace and non-special characters
            pub(super) fn peek_word(&mut self) -> &'a str {
                self.peek_word_internal(false)
            }

            /// First, advances the start cursor past any whitespace,
            /// Then, returns the next contiguos non-whitespace characters,
            /// advancing the start cursor to the end of that returned value
            pub(super) fn eat_word(&mut self) -> &'a str {
                self.peek_word();
                self.take()
            }

            /// First advances the start cursor past any whitespace.
            /// Then, returns a peek of contiguos non-whitespace characters
            pub(super) fn peek_one(&mut self) -> &'a str {
                self.take_whitespace();

                let remaining = &self.data[self.start_cursor..];
                let mut chars = remaining.char_indices().peekable();
                match chars.next() {
                    Some(_) => {}
                    None => {
                        self.end_cursor = self.data.len();
                        return "";
                    }
                }
                match chars.peek() {
                    Some((index, _)) => {
                        self.end_cursor = index + self.start_cursor;
                    }
                    None => {
                        self.end_cursor = self.data.len();
                    }
                }
                &self.data[self.start_cursor..self.end_cursor]
            }

            /// First, advances the start cursor past any whitespace,
            /// Then, returns the next contiguos non-whitespace characters,
            /// advancing the start cursor to the end of that returned value
            pub(super) fn eat_one(&mut self) -> &'a str {
                self.peek_one();
                self.take()
            }

            fn peek_word_internal(&mut self, ignore_special: bool) -> &'a str {
                self.take_whitespace();

                let is_word_delimiter = |c: char| {
                    if ignore_special {
                        Parser::is_word_delimiter(c)
                    } else {
                        Parser::is_word_delimiter(c)
                    }
                };

                let start = self.start_cursor;
                let remaining = &self.data[start..];
                let mut chars = remaining.char_indices().peekable();
                loop {
                    let (_, char) = match chars.next() {
                        Some(c) => c,
                        None => {
                            self.end_cursor = self.data.len();
                            return remaining;
                        }
                    };
                    if !char.is_whitespace() {
                        let (index, is_next_whitespace) = match chars.peek() {
                            Some(c) => (c.0 + start, is_word_delimiter(c.1)),
                            None => (self.data.len(), true),
                        };
                        if is_next_whitespace {
                            self.end_cursor = index;
                            return &self.data[self.start_cursor..self.end_cursor];
                        }
                    }
                    if char.is_whitespace() {
                        chars.peek().map(|c| self.start_cursor = start + c.0);
                    }
                }
            }
            /// First, advances the start cursor past any whitespace.
            /// Then, if the next character is not a quote, returns an empty string.
            /// Otherwise, returns a peek of the set of contiguous characters up to the next quote (including quotes).
            pub(super) fn peek_quoted_string(&mut self) -> &'a str {
                self.take_whitespace();
                let remaining = &self.data[self.start_cursor..];
                let mut chars = remaining.char_indices().peekable();

                match chars.next() {
                    Some(start_quote) if start_quote.1 == '"' => {}
                    _ => {
                        self.end_cursor = self.start_cursor;
                        return "";
                    }
                }

                loop {
                    match chars.next() {
                        Some((_, c)) if !(c == '"') => {}
                        _ => match chars.peek() {
                            Some((index, _)) => {
                                self.end_cursor = index + self.start_cursor;
                                return &self.data[self.start_cursor..self.end_cursor];
                            }
                            _ => {
                                self.end_cursor = self.data.len();
                                return &self.data[self.start_cursor..self.end_cursor];
                            }
                        },
                    }
                }
            }

            pub(super) fn eat_quoted_string(&mut self) -> &'a str {
                self.peek_quoted_string();
                self.take()
            }
        }

        mod test {
            #[test]
            fn test_word() {
                let data = "word";
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_word(), "word");
                assert_eq!(parser.peek_word(), "word");
                assert_eq!(parser.eat_word(), "word");

                let data = "   word";
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_word(), "word");
                assert_eq!(parser.peek_word(), "word");
                assert_eq!(parser.eat_word(), "word");

                let data = "  word   two  ";
                let mut parser: crate::bed::autosql::parse::parser::Parser<'_> =
                    super::Parser::of(data);
                assert_eq!(parser.peek_word(), "word");
                assert_eq!(parser.peek_word(), "word");
                assert_eq!(parser.eat_word(), "word");
                assert_eq!(parser.peek_word(), "two");
                assert_eq!(parser.peek_word(), "two");
                assert_eq!(parser.eat_word(), "two");

                let data = "   word;";
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_word(), "word");
                assert_eq!(parser.peek_word(), "word");
                assert_eq!(parser.eat_word(), "word");
            }

            #[test]
            fn test_quoted_string() {
                let data = r#""simple""#;
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_quoted_string(), r#""simple""#);
                assert_eq!(parser.peek_quoted_string(), r#""simple""#);
                assert_eq!(parser.eat_quoted_string(), r#""simple""#);

                let data = r#"  "simple""#;
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_quoted_string(), r#""simple""#);
                assert_eq!(parser.peek_quoted_string(), r#""simple""#);
                assert_eq!(parser.eat_quoted_string(), r#""simple""#);

                let data = r#"  "simple"  "#;
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_quoted_string(), r#""simple""#);
                assert_eq!(parser.peek_quoted_string(), r#""simple""#);
                assert_eq!(parser.eat_quoted_string(), r#""simple""#);

                let data = r#"  "simple"#;
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_quoted_string(), r#""simple"#);
                assert_eq!(parser.peek_quoted_string(), r#""simple"#);
                assert_eq!(parser.eat_quoted_string(), r#""simple"#);

                let data = r#"  simple"#;
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_quoted_string(), r#""#);
                assert_eq!(parser.peek_quoted_string(), r#""#);
                assert_eq!(parser.eat_quoted_string(), r#""#);

                let data = r#"  "simple"word"#;
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_quoted_string(), r#""simple""#);
                assert_eq!(parser.peek_quoted_string(), r#""simple""#);
                assert_eq!(parser.eat_quoted_string(), r#""simple""#);

                let data = r#""multiple words""#;
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_quoted_string(), r#""multiple words""#);
                assert_eq!(parser.peek_quoted_string(), r#""multiple words""#);
                assert_eq!(parser.eat_quoted_string(), r#""multiple words""#);

                let data = r#"  "multiple words"   test"#;
                let mut parser = super::Parser::of(data);
                assert_eq!(parser.peek_quoted_string(), r#""multiple words""#);
                assert_eq!(parser.peek_quoted_string(), r#""multiple words""#);
                assert_eq!(parser.eat_quoted_string(), r#""multiple words""#);
            }
        }
    }

    #[derive(Debug)]
    pub enum ParseError {
        InvalidDeclareType(String),
        InvalidDeclareName(String),
        InvalidDeclareBrackets(String),
        InvalidFieldSizeClose(String),
        InvalidFieldCommentSeparater(String),
        InvalidFieldValuesBrackets(String),
        InvalidIndexSizeBrackets(String),
    }

    #[derive(Copy, Clone, Debug)]
    pub enum DeclarationType {
        Simple,
        Object,
        Table,
    }

    #[derive(Clone, Debug)]
    pub enum IndexType {
        Primary,
        Index(Option<String>),
        Unique,
    }

    #[derive(Clone, Debug)]
    pub struct DeclareName {
        pub name: String,
        pub index_type: Option<IndexType>,
        pub auto: bool,
    }

    #[derive(Clone, Debug)]
    pub struct Declaration {
        pub declaration_type: DeclarationType,
        pub name: DeclareName,
        pub comment: String,
        pub fields: Vec<Field>,
    }

    #[derive(Clone, Debug)]
    pub enum FieldType {
        Int,
        Uint,
        Short,
        Ushort,
        Byte,
        Ubyte,
        Float,
        Double,
        Char,
        String,
        Lstring,
        Bigint,
        Enum(Vec<String>),
        Set(Vec<String>),
        Declaration(DeclarationType, DeclareName),
    }

    impl FieldType {
        fn try_parse(parser: &mut parser::Parser<'_>) -> Result<Option<Self>, ParseError> {
            let field_type: &str = &parser.peek_word().to_lowercase();
            let field_type = match field_type {
                "int" => FieldType::Int,
                "uint" => FieldType::Uint,
                "short" => FieldType::Short,
                "ushort" => FieldType::Ushort,
                "byte" => FieldType::Byte,
                "ubyte" => FieldType::Ubyte,
                "float" => FieldType::Float,
                "double" => FieldType::Double,
                "char" => FieldType::Char,
                "string" => FieldType::String,
                "lstring" => FieldType::Lstring,
                "bigint" => FieldType::Bigint,
                "enum" => {
                    parser.take();
                    let open_bracket = parser.eat_one();
                    if open_bracket != "(" {
                        return Err(ParseError::InvalidFieldValuesBrackets(
                            open_bracket.to_string(),
                        ));
                    }
                    let mut values = vec![];
                    loop {
                        let value = parser.eat_word();
                        if value == ")" {
                            break;
                        }
                        values.push(value.to_string());
                        let close = parser.eat_one();
                        if close == ")" {
                            break;
                        }
                    }
                    return Ok(Some(FieldType::Enum(values)));
                }
                "set" => {
                    parser.take();
                    let open_bracket = parser.eat_one();
                    if open_bracket != "(" {
                        return Err(ParseError::InvalidFieldValuesBrackets(
                            open_bracket.to_string(),
                        ));
                    }
                    let mut values = vec![];
                    loop {
                        let value = parser.eat_word();
                        if value == ")" {
                            break;
                        }
                        values.push(value.to_string());
                        let close = parser.eat_one();
                        if close == ")" {
                            break;
                        }
                    }
                    return Ok(Some(FieldType::Set(values)));
                }
                "simple" => {
                    parser.take();
                    let declare_name = parse_declare_name(parser)?;
                    return Ok(Some(FieldType::Declaration(
                        DeclarationType::Simple,
                        declare_name,
                    )));
                }
                "object" => {
                    parser.take();
                    let declare_name = parse_declare_name(parser)?;
                    return Ok(Some(FieldType::Declaration(
                        DeclarationType::Object,
                        declare_name,
                    )));
                }
                "table" => {
                    parser.take();
                    let declare_name = parse_declare_name(parser)?;
                    return Ok(Some(FieldType::Declaration(
                        DeclarationType::Object,
                        declare_name,
                    )));
                }
                _ => return Ok(None),
            };
            parser.take();
            return Ok(Some(field_type));
        }
    }

    #[derive(Clone, Debug)]
    pub struct Field {
        pub field_type: FieldType,
        pub field_size: Option<String>,
        pub name: String,
        pub index_type: Option<IndexType>,
        pub auto: bool,
        pub comment: String,
    }

    pub fn parse_autosql(data: &str) -> Result<Vec<Declaration>, ParseError> {
        let mut parser = parser::Parser::of(data);

        parse_declaration_list(&mut parser)
    }

    fn parse_declaration_list(
        parser: &mut parser::Parser<'_>,
    ) -> Result<Vec<Declaration>, ParseError> {
        let mut declarations = vec![];

        let mut i = 0;
        loop {
            if i > 3 {
                break;
            }
            i += 1;
            let dec = parse_declaration(parser)?;
            match dec {
                Some(d) => declarations.push(d),
                None => break,
            }
        }

        Ok(declarations)
    }

    fn parse_declaration(
        parser: &mut parser::Parser<'_>,
    ) -> Result<Option<Declaration>, ParseError> {
        let declare_type = parser.eat_word();
        let declaration_type = match declare_type {
            "simple" => DeclarationType::Simple,
            "object" => DeclarationType::Object,
            "table" => DeclarationType::Table,
            "" => return Ok(None),
            _ => return Err(ParseError::InvalidDeclareType(declare_type.to_string())),
        };

        let declare_name = parse_declare_name(parser)?;

        let comment = parser.eat_quoted_string().to_string();

        let opening_bracket = parser.eat_one();

        if opening_bracket != "(" {
            return Err(ParseError::InvalidDeclareBrackets(
                opening_bracket.to_string(),
            ));
        }

        let fields = parse_field_list(parser)?;

        let closing_bracket = parser.eat_one();

        if closing_bracket != ")" {
            return Err(ParseError::InvalidDeclareBrackets(
                closing_bracket.to_string(),
            ));
        }

        Ok(Some(Declaration {
            declaration_type,
            name: declare_name,
            comment,
            fields,
        }))
    }

    fn parse_declare_name(parser: &mut parser::Parser<'_>) -> Result<DeclareName, ParseError> {
        let declare_name = parser.eat_word();
        if !declare_name.chars().next().unwrap_or(' ').is_alphabetic()
            || declare_name.chars().any(|c| !c.is_alphanumeric())
        {
            return Err(ParseError::InvalidDeclareName(declare_name.to_string()));
        }
        let declare_name = declare_name.to_string();

        let next_word = parser.peek_word();
        let index_type = match next_word {
            "primary" => {
                parser.eat_word();
                Some(IndexType::Primary)
            }
            "index" => {
                parser.eat_word();

                let next = parser.peek_one();
                let size = if next == "[" {
                    parser.eat_one();
                    let size = parser.eat_word().to_string();
                    let close = parser.eat_one();
                    if close != "]" {
                        return Err(ParseError::InvalidIndexSizeBrackets(close.to_string()));
                    }
                    Some(size)
                } else {
                    None
                };
                Some(IndexType::Index(size))
            }
            "unique" => {
                parser.eat_word();
                Some(IndexType::Unique)
            }
            "auto" => None,
            _ => None,
        };

        let next_word = parser.peek_word();
        let auto = if next_word == "auto" {
            parser.eat_word();
            true
        } else {
            false
        };
        Ok(DeclareName {
            name: declare_name,
            index_type,
            auto,
        })
    }

    fn parse_field_list(parser: &mut parser::Parser<'_>) -> Result<Vec<Field>, ParseError> {
        let mut fields = vec![];
        loop {
            let field_type = match FieldType::try_parse(parser)? {
                Some(field_type) => field_type,
                None => break,
            };

            let next_word = parser.peek_one();

            let (field_size, field_name) = if next_word == "[" {
                parser.eat_one();
                let size = parser.eat_word();
                let close = parser.eat_one();
                if close != "]" {
                    return Err(ParseError::InvalidFieldSizeClose(close.to_string()));
                }
                (Some(size.to_string()), parser.eat_word())
            } else {
                let next_word = parser.eat_word();
                (None, next_word)
            };
            let field_name = field_name.to_string();

            let next_word = parser.peek_word();
            let index_type = match next_word {
                "primary" => {
                    parser.eat_word();
                    Some(IndexType::Primary)
                }
                "index" => {
                    parser.eat_word();

                    let next = parser.peek_one();
                    let size = if next == "[" {
                        parser.eat_one();
                        let size = parser.eat_word().to_string();
                        let close = parser.eat_one();
                        if close != "]" {
                            return Err(ParseError::InvalidIndexSizeBrackets(close.to_string()));
                        }
                        Some(size)
                    } else {
                        None
                    };
                    Some(IndexType::Index(size))
                }
                "unique" => {
                    parser.eat_word();
                    Some(IndexType::Unique)
                }
                "auto" => None,
                _ => None,
            };

            let next_word = parser.peek_word();
            let auto = if next_word == "auto" {
                parser.eat_word();
                true
            } else {
                false
            };

            let semicolon = parser.eat_one();
            if semicolon != ";" {
                return Err(ParseError::InvalidFieldCommentSeparater(
                    semicolon.to_string(),
                ));
            }

            let comment = parser.eat_quoted_string().to_string();

            fields.push(Field {
                field_type,
                field_size,
                name: field_name,
                index_type,
                auto,
                comment,
            });

            if parser.peek_one() == ")" {
                break;
            }
        }
        return Ok(fields);
    }

    mod test {
        #[test]
        fn test_bed3() {
            super::parse_autosql(super::super::BED3).unwrap();
        }

        #[test]
        fn test_maintest() {
            let main_test = r#"
            simple pt
            "Two dimensional point"
                (
                int x; "x coor"
                int y; "y coor"
                )
            
            object point
            "Three dimensional point"
                (
                char[12] acc;    "GenBank Accession sequence"
                int x;  "x coor"
                int y;  "y coor"
                int z;  "z coor"
                simple pt pt;  "Transformed point."
                )
            
            table polygon
            "A face"
                (
                uint id;    "Unique ID"
                int pointCount;         "point count"
                object point[pointCount] points;	"Points list"
                simple pt[pointCount] persp;  "Points after perspective transformation"
                )
            
            table polyhedron
            "A 3-d object"
                (
                uint id;   "Unique ID"
                string[2] names;   "Name of this figure"
                int polygonCount;   "Polygon count"
                object polygon[polygonCount] polygons; "Polygons"
                simple pt[2] screenBox; "Bounding box in screen coordinates"
                )
            
            simple twoPoint
            "Two points back to back"
                (
                char[12] name; "name of points"
                simple pt a;  "point a"
                simple pt b; "point b"
                simple pt[2] points; "points as array"
                )
            
            table stringArray
            "An array of strings"
                (
                short numNames;	"Number of names"
                string[numNames] names;   "Array of names"
                )
            "#;

            super::parse_autosql(main_test).unwrap();
        }
        #[test]
        fn test_hardtest() {
            let hard_test = r#"
            object point
            "Three dimensional point"
                (
                int x;  "x coor"
                int y;  "y coor"
                int z;  "z coor"
                )
            
            table autoTest
            "Just a test table"
                (
                uint id; "Unique ID"
                char[12] shortName; "12 character or less name"
                string longName; "full name"
                string[3] aliases; "three nick-names"
                object point threeD;  "Three dimensional coordinate"
                int ptCount;  "number of points"
                short[ptCount] pts;  "point list"
                int difCount;  "number of difs"
                Ubyte [difCount] difs; "dif list"
                int[2] xy;  "2d coordinate"
                int valCount; "value count"
                string[valCount] vals; "list of values"
                double dblVal; "double value"
                float fltVal; "float value"
                double[valCount] dblArray; "double array"
                float[valCount] fltArray; "float array"
                )
            "#;

            super::parse_autosql(hard_test).unwrap();
        }

        #[test]
        fn test_indextest() {
            let index_test = r#"
            table addressBook
            "A simple address book"
                (
                uint id primary auto; "Auto increment primary key"
                string name unique;  "Name - first or last or both, we don't care"
                string address;  "Street address"
                string city index[12];  "City"
                uint zipCode index;  "A zip code is always positive, so can be unsigned"
                char[2] state index;  "Just store the abbreviation for the state"
                )
            
            table symbolCols
            "example of enum and set symbolic columns"
                (
                int id primary auto;                                          "unique id"
                enum(male, female) sex index;                          "enumerated column"
                set(cProg,javaProg,pythonProg,awkProg) skills;   "set column"
                )
            "#;

            super::parse_autosql(index_test).unwrap();
        }
    }
}

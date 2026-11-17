#!/usr/bin/env python3
import os
import re

CHAPTER_MAP = [
    ("source_md_backup/index.md", "index.qmd"),
    ("source_md_backup/00-author.md", "00-author.qmd"),
    ("source_md_backup/01-intro2Genomics.md", "01-intro2Genomics.qmd"),
    ("source_md_backup/02-intro2R.md", "02-intro2R.qmd"),
    ("source_md_backup/03-statsForGenomics.md", "03-statsForGenomics.qmd"),
    ("source_md_backup/04-unsupervisedLearning.md", "04-unsupervisedLearning.qmd"),
    ("source_md_backup/05-supervisedLearning.md", "05-supervisedLearning.qmd"),
    ("source_md_backup/06-genomicIntervals.md", "06-genomicIntervals.qmd"),
    ("source_md_backup/07-Read_Processing.md", "07-Read_Processing.qmd"),
    ("source_md_backup/08-rna-seq-analysis.md", "08-rna-seq-analysis.qmd"),
    ("source_md_backup/09-chip-seq-analysis.md", "09-chip-seq-analysis.qmd"),
    ("source_md_backup/10-bs-seq-analysis.md", "10-bs-seq-analysis.qmd"),
    ("source_md_backup/11-multiomics-analysis.md", "11-multiomics-analysis.qmd"),
    ("source_md_backup/12-references.md", "12-references.qmd"),
]

SEC_IDS = [
    "intro", "Rintro", "stats", "unsupervisedLearning", "supervisedLearning",
    "genomicIntervals", "processingReads", "rnaseqanalysis", "chipseq", "bsseq", "multiomics"
]

def convert_callouts(text):
    for callout in ["tip", "note", "warning", "important", "caution"]:
        pattern = rf'\\BeginKnitrBlock{{rmd{callout}}}<div class="rmd{callout}">(.*?)</div>\\EndKnitrBlock{{rmd{callout}}}'
        def repl(m):
            inner = m.group(1).strip()
            return f"\n\n::: {{.callout-{callout}}}\n{inner}\n:::\n\n"
        text = re.sub(pattern, repl, text, flags=re.DOTALL)
    
    pattern_fr = r'\\BeginKnitrBlock{flushright}<p class="flushright">(.*?)</p>\\EndKnitrBlock{flushright}'
    text = re.sub(pattern_fr, r'<div style="text-align: right;">\n\1\n</div>', text, flags=re.DOTALL)
    return text

def convert_figures(text):
    ref_facet = "Combining two plots using `ggplot2::facet_grid()` function."
    ref_cow = "Combining a histogram and scatter plot using `cowplot` package. The plots are labeled as A and B using the arguments in `plot_grid()` function."
    ref_pwm = "PWM scanning principle. A genomic sequence is scanned by a PWM matrix. This matrix is used to measure how likely it is that the transcription factor will bind each nucleotide in each position. Here we are looking at how likely it is that our TF will bind to the sequence ACACT. The score for this sequence is -3.6. The maximal score obtainable by the PWM is 7.2 and minimum is -5.6. After min-max rescaling, -3.6 corresponds to a 15% score and ACACT is not considered a hit."
    
    text = text.replace("(ref:reffacetHistChp3)", ref_facet)
    text = text.replace("(ref:refcowPlotChp3)", ref_cow)
    text = text.replace("(ref:FigurePWMScanning)", ref_pwm)
    
    fig_pattern = r'<div class="figure"[^>]*>\s*<img\s+src=[\"\x27]([^\"\x27]+)[\"\x27](?:\s+alt=[\"\x27]([^\"\x27]*)[\"\x27])?(?:\s+width=[\"\x27]([^\"\x27]+)[\"\x27])?[^>]*/>\s*(?:<p\s+class=[\"\x27]caption[\"\x27]>(.*?)</p>)?\s*</div>'
    
    def fig_repl(m):
        src = m.group(1)
        alt = m.group(2) or ""
        width = m.group(3) or ""
        cap_raw = m.group(4) or ""
        
        fig_id = ""
        caption = ""
        
        if cap_raw:
            id_m = re.search(r'\(\\\#fig:([a-zA-Z0-9_ -]+)\)', cap_raw)
            if id_m:
                fig_id = id_m.group(1).replace(" ", "-")
                caption = re.sub(r'\(\\\#fig:[a-zA-Z0-9_ -]+\)', '', cap_raw).strip()
            else:
                caption = cap_raw.strip()
        elif alt:
            caption = alt.strip()
            
        attr_parts = []
        if fig_id:
            attr_parts.append(f"#fig-{fig_id}")
        if width:
            attr_parts.append(f'width="{width}"')
        attr_parts.append('fig-align="center"')
        
        attrs = "{" + " ".join(attr_parts) + "}"
        return f"\n\n![{caption}]({src}){attrs}\n\n"
        
    text = re.sub(fig_pattern, fig_repl, text, flags=re.DOTALL)
    return text

def convert_equations(text):
    pattern = r'\\begin\{equation\}\s*(.*?)\s*\(\\\#eq:([a-zA-Z0-9_-]+)\)\s*\\end\{equation\}'
    def eq_repl(m):
        formula = m.group(1).strip()
        eq_id = m.group(2)
        return f"\n$$\n{formula}\n$$ {{#eq-{eq_id}}}\n"
    text = re.sub(pattern, eq_repl, text, flags=re.DOTALL)
    return text

def clean_misc(text):
    text = re.sub(r'\\index\{[^}]*\}', '', text)
    text = re.sub(r'\\frontmatter', '', text)
    text = re.sub(r'^(#+\s+.*?)(\s+\{-\})', r'\1 {.unnumbered}', text, flags=re.MULTILINE)
    return text

def convert_references(text):
    # Fix comma citations: [@citation1, @citation2] -> [@citation1; @citation2]
    text = re.sub(r'\[(@[a-zA-Z0-9_:-]+),\s*(@[a-zA-Z0-9_:-]+)\]', r'[\1; \2]', text)
    
    # Chapter crossrefs:
    for sec_id in SEC_IDS:
        text = re.sub(rf'Chapters?\s+\\@ref\({sec_id}\)', f'@sec-{sec_id}', text)
        text = re.sub(rf'chapters?\s+\\@ref\({sec_id}\)', f'@sec-{sec_id}', text)
        text = re.sub(rf'\\@ref\({sec_id}\)', f'@sec-{sec_id}', text)

    # Figure refs:
    text = re.sub(r'Figures?\s+\\@ref\(fig:([a-zA-Z0-9_ -]+)\)\s+and\s+\\@ref\(fig:([a-zA-Z0-9_ -]+)\)', 
                  lambda m: f"@fig-{m.group(1).replace(' ', '-')} and @fig-{m.group(2).replace(' ', '-')}", text)
    text = re.sub(r'figures?\s+\\@ref\(fig:([a-zA-Z0-9_ -]+)\)\s+and\s+\\@ref\(fig:([a-zA-Z0-9_ -]+)\)', 
                  lambda m: f"@fig-{m.group(1).replace(' ', '-')} and @fig-{m.group(2).replace(' ', '-')}", text)
    text = re.sub(r'Figures?\s+\\@ref\(fig:([a-zA-Z0-9_ -]+)\)', lambda m: f"@fig-{m.group(1).replace(' ', '-')}", text)
    text = re.sub(r'figures?\s+\\@ref\(fig:([a-zA-Z0-9_ -]+)\)', lambda m: f"@fig-{m.group(1).replace(' ', '-')}", text)
    text = re.sub(r'\(\\@ref\(fig:([a-zA-Z0-9_ -]+)\)\)', lambda m: f"(@fig-{m.group(1).replace(' ', '-')})", text)
    text = re.sub(r'\\@ref\(fig:([a-zA-Z0-9_ -]+)\)', lambda m: f"@fig-{m.group(1).replace(' ', '-')}", text)
    
    # Equation refs:
    text = re.sub(r'Equations?\s+\\@ref\(eq:([a-zA-Z0-9_-]+)\)', r'@eq-\1', text)
    text = re.sub(r'Eq\.\s+\\@ref\(eq:([a-zA-Z0-9_-]+)\)', r'@eq-\1', text)
    text = re.sub(r'\\@ref\(eq:([a-zA-Z0-9_-]+)\)', r'@eq-\1', text)
    
    # Cleanup any leftover raw \@ref(ID)
    text = re.sub(r'\\@ref\(([a-zA-Z0-9_-]+)\)', r'@sec-\1', text)
    return text

def process_file(src, dst):
    with open(src, "r", encoding="utf-8") as f:
        content = f.read()

    # Set sec- ID on main chapter title
    for sec_id in SEC_IDS:
        content = re.sub(rf'^(#\s+.*?)\{{#{sec_id}\}}', rf'\1{{#sec-{sec_id}}}', content, flags=re.MULTILINE)

    if "index.md" in src:
        if content.startswith("---"):
            parts = content.split("---", 2)
            if len(parts) >= 3:
                content = parts[2].lstrip()
        if not content.startswith("# Preface"):
            content = "# Preface {.unnumbered}\n\n" + content
        content = content.replace("websites , video", "websites, video")
        content = content.replace("Janathan Kitt", "Jonathan Kitt")

    elif "00-author.md" in src:
        content = content.replace("# About the Authors {-}", "# About the Authors {.unnumbered}")
        content = content.replace("# About the Author(s) {-}", "# About the Authors {.unnumbered}")

    elif "01-intro2Genomics.md" in src:
        content = content.replace("applications in the clinic are also gaining momemntum:", 
                                  "applications in the clinic are also gaining momentum:")
        # Replace histoneMod table with clean markdown table and caption below
        old_histone = """Modifications                 Effect   
------------                  ------  
    H3K9ac                    Active promoters and enhancers        
    H3K14ac                   Active transcription       
    H3K4me3/me2/me1           Active promoters and enhancers, 
                              H3K4me1 and H3K27ac is enhancer-specific
    H3K27ac                   H3K27ac is enhancer-specific          
    H3K36me3                  Active transcribed regions
    H3K27me3/me2/me1          Silent promoters
    H3K9me3/me2/me1           Silent promoters

Table: (\\#tab:histoneMod) Histone modifications and their effects. If more than one histone modification has the same effect, they are separated by commas."""
        new_histone = """| Modifications | Effect |
|:--------------|:-------|
| H3K9ac | Active promoters and enhancers |
| H3K14ac | Active transcription |
| H3K4me3/me2/me1 | Active promoters and enhancers, H3K4me1 and H3K27ac is enhancer-specific |
| H3K27ac | H3K27ac is enhancer-specific |
| H3K36me3 | Active transcribed regions |
| H3K27me3/me2/me1 | Silent promoters |
| H3K9me3/me2/me1 | Silent promoters |

: Histone modifications and their effects. If more than one histone modification has the same effect, they are separated by commas. {#tbl-histonemod}"""
        content = content.replace(old_histone, new_histone)
        content = re.sub(r'Tables?\s+\\@ref\(tab:histoneMod\)', '@tbl-histonemod', content)

    elif "02-intro2R.md" in src:
        content = content.replace("df.f=d(enhancerFilePath, header = FALSE,data.table=FALSE)",
                                  "df.f=fread(enhancerFilePath, header = FALSE,data.table=FALSE)")
        content = content.replace("dim(m1)              # 2 by 5 matrix",
                                  "dim(m1)              # 4 by 2 matrix")
        content = content.replace("interquartile tange", "interquartile range")
        content = content.replace("save one object at a type", "save one object at a time")
        content = content.replace("Scatter plots are one  the most common", "Scatter plots are one of the most common")
        content = content.replace("and plotted those against", "and plot those against")
        content = content.replace("examples of those that in following chapters", "examples of those in following chapters")
        content = content.replace("y itself has not been unchanged", "y itself is unchanged")

    elif "03-statsForGenomics.md" in src:
        content = re.sub(r'the `median\(\)` function\s+will calculate the mean of a provided vector of numbers\.\s*\\index\{median\}',
                         'the `median()` function will calculate the median of a provided vector of numbers.', content)
        content = content.replace("will calculate the mean of a provided vector of numbers.",
                                  "will calculate the median of a provided vector of numbers.")
        content = content.replace("what we want is the are under the curve",
                                  "what we want is the area under the curve")
        content = content.replace("use software for that tables are outdated",
                                  "use software for that since tables are outdated")
        content = content.replace("data=matrix(rnorm(6000,mean=200,sd=70),ncol=6)",
                                  "data=matrix(rnorm(6000,mean=200,sd=70),ncol=6) keeping the number of samples per dataset constant, as n=6")

    elif "04-unsupervisedLearning.md" in src:
        old_exercise = "2. Our next tasks are to remove eigenvectors and reconstruct the matrix using SVD, then calculate the reconstruction error as the difference between original and reconstructed matrix. HINT: You have to use the `svd()` function and equalize eigenvalue to $0$ for the component you want to remove. [Difficulty: **Intermediate/Advanced**]"
        new_exercise = "2. Our next tasks are to remove eigenvectors and reconstruct the matrix using SVD, then calculate the reconstruction error as the difference between original and reconstructed matrix. Remove a few eigenvectors, reconstruct the matrix and calculate the reconstruction error. Reconstruction error can be euclidean distance between original and reconstructed matrices. HINT: You have to use the `svd()` function and equalize eigenvalue to $0$ for the component you want to remove. [Difficulty: **Intermediate/Advanced**]"
        content = content.replace(old_exercise, new_exercise)
        # Table expTable
        old_exptbl = """Table: (\\#tab:expTable)Gene expressions from patients

            IRX4   OCT4   PAX6
---------  -----  -----  -----
patient1      11     10      1
patient2      13     13      3
patient3       2      4     10
patient4       1      3      9"""
        new_exptbl = """|          | IRX4 | OCT4 | PAX6 |
|:---------|-----:|-----:|-----:|
| patient1 |   11 |   10 |    1 |
| patient2 |   13 |   13 |    3 |
| patient3 |    2 |    4 |   10 |
| patient4 |    1 |    3 |    9 |

: Gene expressions from patients {#tbl-exptable}"""
        content = content.replace(old_exptbl, new_exptbl)
        content = re.sub(r'Tables?\s+\\@ref\(tab:expTable\)', '@tbl-exptable', content)

    elif "05-supervisedLearning.md" in src:
        content = content.replace("where $j$ is the number of parameters/coefficients in the model and $L_{log}$ is the log loss function in Eq.",
                                  "where $p$ is the total number of parameters/coefficients in the model and $L_{log}$ is the log loss function in Eq.")

    elif "07-Read_Processing.md" in src:
        content = content.replace("how many reads overlapping with your promoter set",
                                  "how many reads overlap with your promoter set")
        content = content.replace("Dealing with mapped reads are described previously in chapter",
                                  "Dealing with mapped reads is described in chapter")
        content = content.replace("is a fundemental step", "is a fundamental step")

    elif "08-rna-seq-analysis.md" in src:
        # corrplot2 table
        old_corr = """Table: (\\#tab:corrplot2)Correlation scores between samples

             CASE_1      CASE_2      CTRL_1      CTRL_2
-------  ----------  ----------  ----------  ----------
CASE_1    1.0000000   0.9924606   0.9594011   0.9635760
CASE_2    0.9924606   1.0000000   0.9725646   0.9793835
CTRL_1    0.9594011   0.9725646   1.0000000   0.9879862
CTRL_2    0.9635760   0.9793835   0.9879862   1.0000000"""
        new_corr = """|        | CASE_1 | CASE_2 | CTRL_1 | CTRL_2 |
|:-------|-------:|-------:|-------:|-------:|
| CASE_1 | 1.0000000 | 0.9924606 | 0.9594011 | 0.9635760 |
| CASE_2 | 0.9924606 | 1.0000000 | 0.9725646 | 0.9793835 |
| CTRL_1 | 0.9594011 | 0.9725646 | 1.0000000 | 0.9879862 |
| CTRL_2 | 0.9635760 | 0.9793835 | 0.9879862 | 1.0000000 |

: Correlation scores between samples {#tbl-corrplot2}"""
        content = content.replace(old_corr, new_corr)
        content = re.sub(r'\\@ref\(tab:corrplot2\)', '@tbl-corrplot2', content)

        # GO table
        old_go = """Table: (\\#tab:GOanalysistable)Top GO terms sorted by p-value.

      p.value   term.size   precision  domain   term.name                            
---  --------  ----------  ----------  -------  -------------------------------------
64          0        2740       0.223  CC       plasma membrane part                 
23          0        1609       0.136  BP       ion transport                        
16          0        3656       0.258  BP       regulation of biological quality     
30          0         385       0.042  BP       extracellular structure organization 
34          0        7414       0.452  BP       multicellular organismal process     
78          0        1069       0.090  MF       transmembrane transporter activity   
47          0        1073       0.090  BP       organic acid metabolic process       
5           0         975       0.083  BP       response to drug                     
18          0        1351       0.107  BP       biological adhesion                  
31          0        4760       0.302  BP       system development                   """
        new_go = """| | p.value | term.size | precision | domain | term.name |
|:---|--------:|----------:|----------:|:-------|:-------------------------------------|
| 64 | 0 | 2740 | 0.223 | CC | plasma membrane part |
| 23 | 0 | 1609 | 0.136 | BP | ion transport |
| 16 | 0 | 3656 | 0.258 | BP | regulation of biological quality |
| 30 | 0 | 385 | 0.042 | BP | extracellular structure organization |
| 34 | 0 | 7414 | 0.452 | BP | multicellular organismal process |
| 78 | 0 | 1069 | 0.090 | MF | transmembrane transporter activity |
| 47 | 0 | 1073 | 0.090 | BP | organic acid metabolic process |
| 5 | 0 | 975 | 0.083 | BP | response to drug |
| 18 | 0 | 1351 | 0.107 | BP | biological adhesion |
| 31 | 0 | 4760 | 0.302 | BP | system development |

: Top GO terms sorted by p-value. {#tbl-go-analysis}"""
        content = content.replace(old_go, new_go)
        content = re.sub(r'Tables?\s+\\@ref\(tab:GOanalysistable\)', '@tbl-go-analysis', content)

        # gseaPost1
        old_g1 = """Table: (\\#tab:gseaPost1)Up-regulation statistics

               p.geomean   stat.mean    p.val    q.val   set.size     exp1
------------  ----------  ----------  -------  -------  ---------  -------
top_GO_term       0.0000      7.1994   0.0000   0.0000         32   0.0000
random_set        0.5832     -0.2113   0.5832   0.5832         25   0.5832"""
        new_g1 = """| | p.geomean | stat.mean | p.val | q.val | set.size | exp1 |
|:------------|-----------:|----------:|------:|------:|---------:|-----:|
| top_GO_term | 0.0000 | 7.1994 | 0.0000 | 0.0000 | 32 | 0.0000 |
| random_set | 0.5832 | -0.2113 | 0.5832 | 0.5832 | 25 | 0.5832 |

: Up-regulation statistics {#tbl-gsea-post1}"""
        content = content.replace(old_g1, new_g1)
        content = re.sub(r'Tables?\s+\\@ref\(tab:gseaPost1\)', '@tbl-gsea-post1', content)
        content = re.sub(r'\\@ref\(tab:gseaPost1\)', '@tbl-gsea-post1', content)

        # gseaPost2
        old_g2 = """Table: (\\#tab:gseaPost2)Down-regulation statistics

               p.geomean   stat.mean    p.val    q.val   set.size     exp1
------------  ----------  ----------  -------  -------  ---------  -------
random_set        0.4168     -0.2113   0.4168   0.8336         25   0.4168
top_GO_term       1.0000      7.1994   1.0000   1.0000         32   1.0000"""
        new_g2 = """| | p.geomean | stat.mean | p.val | q.val | set.size | exp1 |
|:------------|-----------:|----------:|------:|------:|---------:|-----:|
| random_set | 0.4168 | -0.2113 | 0.4168 | 0.8336 | 25 | 0.4168 |
| top_GO_term | 1.0000 | 7.1994 | 1.0000 | 1.0000 | 32 | 1.0000 |

: Down-regulation statistics {#tbl-gsea-post2}"""
        content = content.replace(old_g2, new_g2)
        content = re.sub(r'Tables?\s+\\@ref\(tab:gseaPost2\)', '@tbl-gsea-post2', content)
        content = re.sub(r'\\@ref\(tab:gseaPost2\)', '@tbl-gsea-post2', content)

    elif "09-chip-seq-analysis.md" in src:
        old_pattern = '2. **nmotifs** - the number of motifs to look for.\n\n\n\n```\n## top 3  4, 5-mers: 12 40 52'
        new_code = '''2. **nmotifs** - the number of motifs to look for.

```r
# run the rGADEM motif discovery
novel_motifs = GADEM(
  Sequences = ctcf_seq,
  seed = 1000,
  nmotifs = 5)
```

```
## top 3  4, 5-mers: 12 40 52'''
        content = content.replace(old_pattern, new_code)

        # peakqualityshow table
        old_ctcf = """Table: (\\#tab:peakqualityshow)Position Frequency Matrix (PFM) for the CTCF motif

         1      2      3      4      5      6      7      8      9   10     11     12   13     14     15     16     17   18     19     20
---  -----  -----  -----  -----  -----  -----  -----  -----  -----  ---  -----  -----  ---  -----  -----  -----  -----  ---  -----  -----
A     0.17   0.23   0.29   0.10   0.33   0.06   0.05   0.04   0.02    0   0.25   0.00    0   0.05   0.25   0.00   0.17    0   0.02   0.19
C     0.42   0.28   0.30   0.32   0.11   0.33   0.56   0.00   0.96    1   0.67   0.69    1   0.04   0.07   0.42   0.15    0   0.06   0.43
G     0.25   0.23   0.26   0.27   0.42   0.55   0.05   0.83   0.01    0   0.03   0.00    0   0.02   0.53   0.55   0.05    1   0.87   0.15
T     0.16   0.27   0.15   0.31   0.14   0.06   0.33   0.13   0.00    0   0.06   0.31    0   0.89   0.15   0.03   0.62    0   0.05   0.23"""
        new_ctcf = """| | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
|:---|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|---:|-----:|-----:|---:|-----:|-----:|-----:|-----:|---:|-----:|-----:|
| A | 0.17 | 0.23 | 0.29 | 0.10 | 0.33 | 0.06 | 0.05 | 0.04 | 0.02 | 0 | 0.25 | 0.00 | 0 | 0.05 | 0.25 | 0.00 | 0.17 | 0 | 0.02 | 0.19 |
| C | 0.42 | 0.28 | 0.30 | 0.32 | 0.11 | 0.33 | 0.56 | 0.00 | 0.96 | 1 | 0.67 | 0.69 | 1 | 0.04 | 0.07 | 0.42 | 0.15 | 0 | 0.06 | 0.43 |
| G | 0.25 | 0.23 | 0.26 | 0.27 | 0.42 | 0.55 | 0.05 | 0.83 | 0.01 | 0 | 0.03 | 0.00 | 0 | 0.02 | 0.53 | 0.55 | 0.05 | 1 | 0.87 | 0.15 |
| T | 0.16 | 0.27 | 0.15 | 0.31 | 0.14 | 0.06 | 0.33 | 0.13 | 0.00 | 0 | 0.06 | 0.31 | 0 | 0.89 | 0.15 | 0.03 | 0.62 | 0 | 0.05 | 0.23 |

: Position Frequency Matrix (PFM) for the CTCF motif {#tbl-peakqualityshow}"""
        content = content.replace(old_ctcf, new_ctcf)
        content = re.sub(r'Tables?\s+\\@ref\(tab:peakqualityshow\)', '@tbl-peakqualityshow', content)

    elif "10-bs-seq-analysis.md" in src:
        content = content.replace("### How DNA methylation is set ?", "### How DNA methylation is set?")
        old_pattern = 'We will first get a list of file paths and have a look at the content.\n\n\n\n\nIf you look what is inside the `file.list` variable'
        new_code = '''We will first get a list of file paths and have a look at the content.

```r
library(methylKit)
file.list = list( system.file("extdata", 
                             "test1.myCpG.txt", package = "methylKit"),
                 system.file("extdata", 
                             "test2.myCpG.txt", package = "methylKit"),
                 system.file("extdata", 
                             "control1.myCpG.txt", package = "methylKit"),
                 system.file("extdata", 
                             "control2.myCpG.txt", package = "methylKit") )
```

If you look what is inside the `file.list` variable'''
        content = content.replace(old_pattern, new_code)

    elif "11-multiomics-analysis.md" in src:
        content = content.replace("A heatmap of NMF factors shows separability of tumors into subtype clusters. This plot is more useful than a scatter plot when there are more than two factors.",
                                  "A heatmap of 2-component NMF factors shows the separability of tumors into subtype clusters. A heatmap representation like this is especially useful when scaling beyond two factors where scatter plots become less practical.")

        # Table GE
        old_ge = """Table: (\\#tab:moloadMultiomicsGE)Example gene expression data (head)

                 RNF113A    S100A13      AP3D1   ATP6V1G1     UBQLN4      TPPP3
-------------  ---------  ---------  ---------  ---------  ---------  ---------
TCGA.A6.2672    21.19567   19.72600   11.53022    0.00000   15.35637   12.76747
TCGA.A6.3809    21.50866   18.65729   12.98830   14.12675   19.62208    0.00000
TCGA.A6.5661    20.08072   18.97034   10.83759   15.31325    0.00000    0.00000
TCGA.A6.5665     0.00000   11.88336   10.24248   19.79300    0.00000    0.00000
TCGA.A6.6653     0.00000   12.07753    0.00000    0.00000    0.00000    0.00000
TCGA.A6.6780     0.00000   12.99128    0.00000   19.96976   13.17618   11.58742"""
        new_ge = """| | RNF113A | S100A13 | AP3D1 | ATP6V1G1 | UBQLN4 | TPPP3 |
|:---|---------:|---------:|---------:|---------:|---------:|---------:|
| TCGA.A6.2672 | 21.19567 | 19.72600 | 11.53022 | 0.00000 | 15.35637 | 12.76747 |
| TCGA.A6.3809 | 21.50866 | 18.65729 | 12.98830 | 14.12675 | 19.62208 | 0.00000 |
| TCGA.A6.5661 | 20.08072 | 18.97034 | 10.83759 | 15.31325 | 0.00000 | 0.00000 |
| TCGA.A6.5665 | 0.00000 | 11.88336 | 10.24248 | 19.79300 | 0.00000 | 0.00000 |
| TCGA.A6.6653 | 0.00000 | 12.07753 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| TCGA.A6.6780 | 0.00000 | 12.99128 | 0.00000 | 19.96976 | 13.17618 | 11.58742 |

: Example gene expression data (head) {#tbl-moload-ge}"""
        content = content.replace(old_ge, new_ge)
        content = re.sub(r'Tables?\s+\\@ref\(tab:moloadMultiomicsGE\)', '@tbl-moload-ge', content)

        # Table MUT
        old_mut = """Table: (\\#tab:moloadMultiomicsMUT)Example mutation data (head)

                TTN   TP53   APC   KRAS   SYNE1   MUC16
-------------  ----  -----  ----  -----  ------  ------
TCGA.A6.2672      1      0     0      0       1       1
TCGA.A6.3809      1      0     0      0       0       0
TCGA.A6.5661      1      0     0      0       1       1
TCGA.A6.5665      1      0     0      0       1       1
TCGA.A6.6653      1      0     0      1       0       0
TCGA.A6.6780      1      0     0      0       0       1"""
        new_mut = """| | TTN | TP53 | APC | KRAS | SYNE1 | MUC16 |
|:---|---:|---:|---:|---:|---:|---:|
| TCGA.A6.2672 | 1 | 0 | 0 | 0 | 1 | 1 |
| TCGA.A6.3809 | 1 | 0 | 0 | 0 | 0 | 0 |
| TCGA.A6.5661 | 1 | 0 | 0 | 0 | 1 | 1 |
| TCGA.A6.5665 | 1 | 0 | 0 | 0 | 1 | 1 |
| TCGA.A6.6653 | 1 | 0 | 0 | 1 | 0 | 0 |
| TCGA.A6.6780 | 1 | 0 | 0 | 0 | 0 | 1 |

: Example mutation data (head) {#tbl-moload-mut}"""
        content = content.replace(old_mut, new_mut)
        content = re.sub(r'Tables?\s+\\@ref\(tab:moloadMultiomicsMUT\)', '@tbl-moload-mut', content)

        # Table CNV
        old_cnv = """Table: (\\#tab:moloadMultiomicsCNV)Example copy number data for CRC samples

                8p23.2   8p23.3   8p23.1   8p21.3   8p12   8p22
-------------  -------  -------  -------  -------  -----  -----
TCGA.A6.2672         0        0        0        0      0      0
TCGA.A6.3809         0        0        0        0      0      0
TCGA.A6.5661         0        0        0        0      0      0
TCGA.A6.5665         0        0        0        0      0      0
TCGA.A6.6653         0        0        0        0      0      0
TCGA.A6.6780         0        0        0        0      0      0"""
        new_cnv = """| | 8p23.2 | 8p23.3 | 8p23.1 | 8p21.3 | 8p12 | 8p22 |
|:---|------:|------:|------:|------:|----:|----:|
| TCGA.A6.2672 | 0 | 0 | 0 | 0 | 0 | 0 |
| TCGA.A6.3809 | 0 | 0 | 0 | 0 | 0 | 0 |
| TCGA.A6.5661 | 0 | 0 | 0 | 0 | 0 | 0 |
| TCGA.A6.5665 | 0 | 0 | 0 | 0 | 0 | 0 |
| TCGA.A6.6653 | 0 | 0 | 0 | 0 | 0 | 0 |
| TCGA.A6.6780 | 0 | 0 | 0 | 0 | 0 | 0 |

: Example copy number data for CRC samples {#tbl-moload-cnv}"""
        content = content.replace(old_cnv, new_cnv)
        content = re.sub(r'tables?\s+\\@ref\(tab:moloadMultiomicsCNV\)', '@tbl-moload-cnv', content)
        content = re.sub(r'\\@ref\(tab:moloadMultiomicsCNV\)', '@tbl-moload-cnv', content)

        # Table split 1
        old_sp1 = """Table: (\\#tab:mocnvsplitcolshow1)Example copy number data. Data can be both positive (amplified regions) or negative (deleted regions).

         seg1   seg2
-------  -----  -----
tumor 1  +2     -1   
tumor 2  -1     0    
tumor 3  0      +2   """
        new_sp1 = """| | seg1 | seg2 |
|:---|:---|:---|
| tumor 1 | +2 | -1 |
| tumor 2 | -1 | 0 |
| tumor 3 | 0 | +2 |

: Example copy number data. Data can be both positive (amplified regions) or negative (deleted regions). {#tbl-mocnv-split1}"""
        content = content.replace(old_sp1, new_sp1)
        content = re.sub(r'Tables?\s+\\@ref\(tab:mocnvsplitcolshow1\)', '@tbl-mocnv-split1', content)
        content = re.sub(r'\\@ref\(tab:mocnvsplitcolshow1\)', '@tbl-mocnv-split1', content)

        # Table split 2
        old_sp2 = """Table: (\\#tab:mocnvsplitcolshow2)Example copy number data after splitting each column into a column representing copy number gains (+) and a column representing deletions (-). This data matrix is non-negative, and thus suitable for NMF algorithms.

         seg1+   seg1-   seg2+   seg2-
-------  ------  ------  ------  ------
tumor 1  2       0       0       1     
tumor 2  0       1       0       0     
tumor 3  0       0       2       0     """
        new_sp2 = """| | seg1+ | seg1- | seg2+ | seg2- |
|:---|:---|:---|:---|:---|
| tumor 1 | 2 | 0 | 0 | 1 |
| tumor 2 | 0 | 1 | 0 | 0 |
| tumor 3 | 0 | 0 | 2 | 0 |

: Example copy number data after splitting each column into a column representing copy number gains (+) and a column representing deletions (-). This data matrix is non-negative, and thus suitable for NMF algorithms. {#tbl-mocnv-split2}"""
        content = content.replace(old_sp2, new_sp2)
        content = re.sub(r'Tables?\s+\\@ref\(tab:mocnvsplitcolshow2\)', '@tbl-mocnv-split2', content)
        content = re.sub(r'\\@ref\(tab:mocnvsplitcolshow2\)', '@tbl-mocnv-split2', content)

        # Table moNMFGOTerms
        latex_tbl = r'\\begin\{table\}\s*\\caption\{\(\\\#tab:moNMFGOTerms\)GO-terms associated with NMF factor 2\}.*?\\end\{table\}'
        md_tbl = """| Term | Adjusted.P.value | Combined.Score |
|:-----|-----------------:|---------------:|
| nuclear-transcribed mRNA catabolic process (GO:0000956) | 0 | 207.7403 |
| rRNA metabolic process (GO:0016072) | 0 | 161.4781 |
| nuclear-transcribed mRNA catabolic process, nonsense-mediated decay (GO:0000184) | 0 | 220.4298 |

: GO-terms associated with NMF factor 2 {#tbl-monmf-goterms}"""
        content = re.sub(latex_tbl, md_tbl, content, flags=re.DOTALL)
        content = re.sub(r'Tables?\s+\\@ref\(tab:moNMFGOTerms\)', '@tbl-monmf-goterms', content)
        content = re.sub(r'\\@ref\(tab:moNMFGOTerms\)', '@tbl-monmf-goterms', content)

    elif "12-references.md" in src:
        content = "# References {.unnumbered}\n\n::: {#refs}\n:::\n"

    # Conversions
    content = convert_callouts(content)
    content = convert_figures(content)
    content = convert_equations(content)
    content = convert_references(content)
    content = clean_misc(content)

    with open(dst, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Converted {src} -> {dst}")

for src, dst in CHAPTER_MAP:
    process_file(src, dst)

print("\nAll chapters converted!")

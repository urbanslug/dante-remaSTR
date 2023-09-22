import gzip
import numpy as np
import re

from src.postfilter import PostFilter

contents = """
<table class="mtg" id="content-tg">
    <thead>
        <tr>
            <th class="mtg-s6z2">Motif</th>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>
"""

make_datatable_string = """
<script>
    $(document).ready( function () {{
    $('#content-tg').DataTable();
}} );
</script>
"""

content_string = """ <tr>
    <td class="mtg-s6z2">
        <button class="tablinks" onclick="openTab(event, '{motif_name}'); $('.{motif_name}').trigger('content-change');">
            <a href="#{motif_name}">{motif}</a>
        </button>
    </td>
</tr>"""

content_string_empty = ""

row_string = """  <tr>
    <td class="tg-s6z2">{motif_name}</td>
    <td class="tg-s6z2">{str_seq}</td>
    <td class="tg-s6z2">{allele1}</td>
    <td class="tg-s6z2">{allele1_conf:.1%}</td>
    <td class="tg-s6z2">{allele2}</td>
    <td class="tg-s6z2">{allele2_conf:.1%}</td>
    <td class="tg-s6z2">{motif_conf:.1%}</td>
    <td class="tg-s6z2">{reads_blue}</td>
    <td class="tg-s6z2">{reads_grey}</td>
    <td class="tg-s6z2">{post_bases}</td>
    <td class="tg-s6z2">{post_reps}</td>
    <td class="tg-s6z2">{post_errors}</td>
    <td class="tg-s6z2">{motif_seq}</td>
  </tr>"""

row_string_empty_old = """  <tr>
    <td class="tg-s6z2">{motif_name}</td>
    <td class="tg-s6z2">{str_seq}</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">{reads_blue}</td>
    <td class="tg-s6z2">{reads_grey}</td>
    <td class="tg-s6z2">{post_bases}</td>
    <td class="tg-s6z2">{post_reps}</td>
    <td class="tg-s6z2">{post_errors}</td>
    <td class="tg-s6z2">{motif_seq}</td>
  </tr>"""

row_string_empty = ""

nomenclature_string = """
<tr>
    <td>{count}</td>
    <td>{ref}</td>
    {parts}
</tr>
"""

motif_summary = """
<div class="tabcontent" id="{motif_id}" style="display: none">
<h2 class="summary_nomenclatures">Nomenclatures</h2>
<table class="nomtg">
    <tbody>
        {nomenclatures}
    </tbody>
</table>
<h2 class="summary">Summary table</h2>
<table class="tg" id="tg-{motif_id}">
    <thead>
        <tr>
            <th class="tg-s6z2" rowspan="2">Motif</th>
            <th class="tg-s6z2" rowspan="2">STR<br>sequence</th>
            <th class="tg-s6z2" colspan="2">Allele 1</th>
            <th class="tg-s6z2" colspan="2">Allele 2</th>
            <th class="tg-s6z2" rowspan="2">Overall<br>confidence</th>
            <th class="tg-s6z2" colspan="2">Reads</th>
            <th class="tg-s6z2" colspan="3">Postfilter</th>
            <th class="tg-s6z2" rowspan="2">Sequence</th>
        </tr>
        <tr>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">full</td>
            <td class="tg-s6z2">partial</td>
            <td class="tg-s6z2">bases</td>
            <td class="tg-s6z2">modules</td>
            <td class="tg-s6z2">max. errors</td>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>

<script>
    $(document).ready( function () {{
    $('#tg-{motif_id}').DataTable();
}} );
</script>
<p><a href="#content">Back to content</a></p>
{motifs}
</div>
"""

motif_summary_static = """
<div style="margin: 25px 0 25px 0">
<h2 class="summary">Summary table</h2>
<table class="tg">
    <thead>
        <tr>
            <th class="tg-s6z2" rowspan="2">Motif</th>
            <th class="tg-s6z2" rowspan="2">STR<br>sequence</th>
            <th class="tg-s6z2" colspan="2">Allele 1</th>
            <th class="tg-s6z2" colspan="2">Allele 2</th>
            <th class="tg-s6z2" rowspan="2">Overall<br>confidence</th>
            <th class="tg-s6z2" colspan="2">Reads</th>
            <th class="tg-s6z2" colspan="3">Postfilter</th>
            <th class="tg-s6z2" rowspan="2">Sequence</th>
        </tr>
        <tr>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">full</td>
            <td class="tg-s6z2">partial</td>
            <td class="tg-s6z2">bases</td>
            <td class="tg-s6z2">modules</td>
            <td class="tg-s6z2">max. errors</td>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>
</div>
"""

motif_string = """<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="2">
            <img class="hist pic100" alt="{motif_name} repetitions" src="{motif_reps}" />
        </td>
        <td colspan="1">
            <img class="pcol pic100" alt="{motif_name} pcolor" src="{motif_pcolor}" />
        </td>
    </tr>
</table>
<p><a href="{alignment}">Link to alignments</a></p>
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64 = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="2">
            <div class="hist pic100 {motif_id}" id="hist-{motif_name}"></div>
            <script>
                {{
                    let hist_data = {motif_reps};
                    $(document).ready( function() {{
                        $('.{motif_id}').bind("content-change", function() {{
                            if (document.getElementById('{motif_id}').style.display === 'block') {{
                                Plotly.react('hist-{motif_name}', hist_data, {{}});
                            }}
                            else {{
                                Plotly.react('hist-{motif_name}', {{}}, {{}});
                            }}
                        }})
                    }})
                }}
            </script>
        </td>
        <td colspan="1">
            <div class="pcol pic100 {motif_id}" id="pcol-{motif_name}"></div>
            <script>
                {{
                    let pcol_data = {motif_pcolor};
                    $(document).ready( function() {{
                        $('.{motif_id}').bind("content-change", function() {{
                            if (document.getElementById('{motif_id}').style.display === 'block') {{
                                Plotly.react('pcol-{motif_name}', pcol_data, {{}});
                            }}
                            else {{
                                Plotly.react('pcol-{motif_name}', {{}}, {{}});
                            }}
                        }})
                    }})
                }}
            </script>
        </td>
    </tr>
</table>
<p><a href="{alignment}">Link to alignments</a></p>
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_reponly = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="1">
            <div class="hist pic50 {motif_id}" id="hist2d-{motif_name}"></div>
            <script>
                {{
                    let hist2d_data = {motif_reps};
                    $(document).ready( function() {{
                        $('.{motif_id}').bind("content-change", function() {{
                            if (document.getElementById('{motif_id}').style.display === 'block') {{
                                Plotly.react('hist2d-{motif_name}', hist2d_data, {{}});
                            }}
                            else {{
                                Plotly.react('hist2d-{motif_name}', {{}}, {{}});
                            }}
                        }})
                    }})
                }}
            </script>
        </td>
    </tr>
</table>
<p><a href="{alignment}">Link to alignments</a></p>
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_static = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="2">
            <div class="hist pic100 {motif_id}" id="hist-{motif_name}"></div>
            <script>
                Plotly.react('hist-{motif_name}', {motif_reps}, {{}});
            </script>
        </td>
        <td colspan="1">
            <div class="pcol pic100 {motif_id}" id="pcol-{motif_name}"></div>
            <script>
                Plotly.react('pcol-{motif_name}', {motif_pcolor}, {{}});
            </script>
        </td>
    </tr>
</table>
<p><a href="{alignment}">Link to alignments</a></p>
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_reponly_static = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="1">
            <div class="hist pic50 {motif_id}" id="hist2d-{motif_name}"></div>
            <script>
                Plotly.react('hist2d-{motif_name}', {motif_reps}, {{}});
            </script>
        </td>
    </tr>
</table>
<p><a href="#content">Back to content</a></p>
"""

motif_string_empty_old = """<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
NOT AVAILABLE
<p><a href="#content">Back to content</a></p>
"""

motif_string_empty = ""

alignment_string = """
  <p>{sequence}</p>
  {alignment}
  <hr>
"""

align_vis = """
  <details>
    <summary>{display_text}</summary>
    <div id="A{name}" class="align">press "Run with JS"</div>
    <script>
        var fasta = `{fasta}`;
        var seqs = msa.io.fasta.parse(fasta);
        var opts = {{
            el: document.getElementById("A{name}"),
            vis: {{
                conserv: false,
                metaIdentity: true,
                overviewbox: true,
                seqlogo: true
            }},
            seqs: seqs,
            colorscheme: {{"scheme": "nucleotide"}},
            // smaller menu for JSBin
            menu: "small",
            bootstrapMenu: true
        }};
        var m = new msa.msa(opts);
        m.render()
    </script>
  </details>
"""


def highlight_subpart(seq: str, highlight: int | list[int]) -> tuple[str, str]:
    """
    Highlights subpart of a motif sequence
    :param seq: str - motif sequence
    :param highlight: int/list(int) - part ot highlight
    :return: str, str - motif sequence with highlighted subpart, highlighted subpart
    """
    if highlight is None:
        return seq, ''

    str_part = []
    highlight = np.array(highlight)
    split = [f'{s}]' for s in seq.split(']') if s != '']
    for h in highlight:
        str_part.append(split[h])
        split[h] = f'<b><u>{split[h]}</u></b>'
    return ''.join(split), ''.join(str_part)


def generate_row(motif: str, sequence: str, confidence: tuple[float, int, int, float, float, float, float, float, float],
                 postfilter: PostFilter, reads_blue: int, reads_grey: int, highlight: list[int] = None):
    """
    Generate rows of a summary table in html report.
    :param motif: str - motif name
    :param sequence: str - motif sequence
    :param confidence: tuple - motif confidences and allele predictions
    :param postfilter: PostFilter - postfilter dict from config
    :param reads_blue: int - number of full reads
    :param reads_grey: int - number of partial reads
    :param highlight: list(int)/None - which part of seq to highlight
    :return: str - html string with rows of the summary table
    """
    sequence, subpart = highlight_subpart(sequence, highlight)

    # shorten sequence:
    keep = 10
    first = sequence.find(',')
    last = sequence.rfind(',')
    smaller_seq = sequence if first == -1 else '...' + sequence[first - keep:last + keep + 1] + '...'

    # errors:
    errors = f'{postfilter.max_rel_error * 100:.0f}%'
    if postfilter.max_abs_error is None:
        errors += f' (abs={postfilter.max_abs_error})'

    # fill templates:
    if confidence is None:
        return row_string_empty_old.format(post_bases=postfilter.min_rep_len, post_reps=postfilter.min_rep_cnt,
                                           motif_name=motif, motif_seq=smaller_seq, reads_blue=reads_blue,
                                           reads_grey=reads_grey, str_seq=subpart, post_errors=errors)
    else:
        (c, a1, a2, c1, c2, _, _, _, _) = confidence
        if a1 == 0 and a2 == 0:
            a1 = 'BG'
            a2 = 'BG'
        return row_string.format(post_bases=postfilter.min_rep_len, post_reps=postfilter.min_rep_cnt, motif_name=motif,
                                 motif_seq=smaller_seq, reads_blue=reads_blue, reads_grey=reads_grey, motif_conf=c,
                                 allele1=a1, allele2=a2, allele1_conf=c1, allele2_conf=c2, str_seq=subpart,
                                 post_errors=errors)


def get_alignment_name(alignment_file: str, allele: int) -> str:
    """
    Get alignment file of subpart of the alignment with allele count specified.
    :param alignment_file: str - alignment file name
    :param allele: int - allele repetition count
    :return: str - alignment file name for the subpart of alignment
    """
    # find where is .fasta
    fasta_index = alignment_file.rfind('.fasta')
    # insert '_aX' before .fasta
    return alignment_file[:fasta_index] + '_a' + str(allele) + alignment_file[fasta_index:]


def generate_motifb64(motif_name: str, description: str, sequence: str, repetition: str, pcolor: str, alignment: str | None,
                      filtered_alignment: str | None, confidence: tuple[float, int, int, float, float, float, float, float, float],
                      postfilter: PostFilter, highlight: list[int] = None, static: bool = False) -> tuple[str, str, tuple[str, str]]:
    """
    Generate part of a html report for each motif.
    :param motif_name: str - motif name
    :param description: str - motif description
    :param sequence: str - motif sequence
    :param repetition: str - filename of repetitions figures
    :param pcolor: str - filename of pcolor figures
    :param alignment: str/None - filename of alignment file
    :param filtered_alignment: str/None - filename of filtered alignment file
    :param confidence: tuple - motif confidences and allele predictions
    :param postfilter: PostFilter - postfilter arguments
    :param highlight: list(int)/None - which part of seq to highlight
    :param static: bool - generate static code?
    :return: (str, str) - content and main part of the html report for motifs
    """
    # prepare and generate alignments
    sequence, subpart = highlight_subpart(sequence, highlight)
    motif = f'{motif_name} &ndash; {description}'
    motif_name = f'{motif_name.replace("/", "_")}_{",".join(map(str, highlight)) if highlight is not None else "mot"}'
    motif_clean = re.sub(r'[^\w_]', '', motif_name)
    motif_clean_id = motif_clean.rsplit('_', 1)[0] if highlight == [1] else motif_clean  # trick to solve static html
    align_html_a1 = ''
    align_html_a2 = ''
    if confidence is None:
        result = '-- (---.-%%) -- (---.-%%) total ---.-%%'
    else:
        (c, a1, a2, c1, c2, _, _, _, _) = confidence
        if (a1 == 'B' and a2 == 'B') or (a1 == 0 and a2 == 0):
            result = f'BG {c * 100: 5.1f}%'
        else:
            result = f'{str(a1):2s} ({c1 * 100: 5.1f}%) {str(a2):2s} ({c2 * 100: 5.1f}%) total {c * 100: 5.1f}%'
            if alignment is not None:
                align_html_a1 = generate_alignment(f'{motif_clean}_{str(a1)}', get_alignment_name(alignment, a1), motif_clean.split('_')[0],
                                                   f'Allele 1 ({str(a1):2s}) alignment visualization')
                if a1 != a2:
                    align_html_a2 = generate_alignment(f'{motif_clean}_{str(a2)}', get_alignment_name(alignment, a2), motif_clean.split('_')[0],
                                                       f'Allele 2 ({str(a2):2s}) alignment visualization')

    # errors:
    errors = f'{postfilter.max_rel_error * 100:.0f}%'
    if postfilter.max_abs_error is None:
        errors += f' (abs={postfilter.max_abs_error})'

    # return content and picture parts:
    motif_templates = {'static': {'pcol': motif_stringb64_static, 'no-pcol': motif_stringb64_reponly_static},
                       'dynamic': {'pcol': motif_stringb64, 'no-pcol': motif_stringb64_reponly}}

    if repetition is not None:
        reps = open(repetition, 'r').read()
        align_html = generate_alignment(motif_clean, alignment, motif_clean.split('_')[0])
        filt_align_html = generate_alignment(motif_clean + '_filtered', filtered_alignment, motif_clean.split('_')[0],
                                             'Partial reads alignment visualization')
        # select template
        motif_template = motif_templates['static' if static else 'dynamic']['no-pcol' if pcolor is None else 'pcol']

        # read pcolor if available
        pcol = '' if pcolor is None else open(pcolor, 'r').read()

        # return filled valid template
        return (content_string.format(motif_name=motif_clean.rsplit('_', 1)[0], motif=motif),
                motif_template.format(post_bases=postfilter.min_rep_len, post_reps=postfilter.min_rep_cnt, motif_name=motif_clean_id,
                                      motif_id=motif_clean.rsplit('_', 1)[0], motif=motif, motif_reps=reps, result=result, motif_pcolor=pcol,
                                      alignment=f'{motif_name}/alignments.html', sequence=sequence, errors=errors),
                (motif, alignment_string.format(sequence=sequence, alignment=align_html + align_html_a1 + align_html_a2 + filt_align_html)))
    else:
        return (content_string_empty.format(motif_name=motif_clean.rsplit('_', 1)[0], motif=motif),
                motif_string_empty.format(post_bases=postfilter.min_rep_len, post_reps=postfilter.min_rep_cnt,
                                          motif_name=motif_clean, motif=motif, sequence=sequence, errors=errors),
                (motif, ''))


def generate_alignment(motif: str, alignment_file: str, motif_id: str, display_text: str = 'Click to toggle alignment visualization'):
    """
    Generate HTML code for the fancy alignment.
    :param motif: str - name of the motif
    :param alignment_file: str - filename of the alignment file
    :param motif_id: str - motif identification
    :param display_text: str - string to display when the alignment is hidden
    :return: str - code of the fancy alignment
    """
    if alignment_file is None:
        return ''

    try:
        with gzip.open(alignment_file, 'rt') if alignment_file.endswith('.gz') else open(alignment_file) as f:
            string = f.read()
        debug = string.find('#')

        return align_vis.format(fasta=string[:debug], name=motif, motif_id=motif_id, display_text=display_text)
    except (IOError, TypeError, AttributeError):
        return ''

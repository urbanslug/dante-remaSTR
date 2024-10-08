from jinja2 import Environment, FileSystemLoader, select_autoescape

motifs: list = [
    {
        "name": "ALS",
        "str": "ALS"
    }
]

#     {
#         "name": "BPEC",
#         "str": "BPEC"
#     },
env = Environment(
    loader=FileSystemLoader(["./dante_remastr_standalone_templates/"]),
    autoescape=select_autoescape()
)
template = env.get_template("report_template.html")

print(template.render(motifs=motifs))

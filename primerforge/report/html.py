"""
HTML report generator.
"""

from pathlib import Path
from datetime import datetime

from primerforge.models.pair import PrimerPair


class HTMLReport:
    """
    Generate an HTML report for PrimerForge results.
    """

    def write(
        self,
        pairs: list[PrimerPair],
        output: Path,
    ):

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        rows = ""

        for i, pair in enumerate(
            pairs,
            start=1,
        ):

            passed = (
                "✅ PASS"
                if pair.passed_specificity
                else "❌ FAIL"
            )

            rows += f"""
<tr>
<td>{i}</td>
<td>{pair.forward.sequence}</td>
<td>{pair.reverse.sequence}</td>
<td>{pair.product_size}</td>
<td>{pair.forward.tm:.2f}</td>
<td>{pair.reverse.tm:.2f}</td>
<td>{pair.forward.gc:.1f}</td>
<td>{pair.reverse.gc:.1f}</td>
<td>{pair.population_conservation:.2f}</td>
<td>{pair.population_coverage:.2f}</td>
<td>{pair.specificity_score:.2f}</td>
<td>{passed}</td>
<td>{pair.score:.2f}</td>
</tr>
"""

        html = f"""
<!DOCTYPE html>
<html>

<head>

<meta charset="utf-8">

<title>PrimerForge Report</title>

<style>

body {{
    font-family: Arial, Helvetica, sans-serif;
    margin:40px;
}}

h1 {{
    color:#1f4e79;
}}

table {{
    border-collapse:collapse;
    width:100%;
}}

th {{
    background:#1f4e79;
    color:white;
    padding:8px;
}}

td {{
    border:1px solid #cccccc;
    padding:8px;
    text-align:center;
}}

tr:nth-child(even){{
    background:#f5f5f5;
}}

.pass {{
    color:green;
    font-weight:bold;
}}

.fail {{
    color:red;
    font-weight:bold;
}}

</style>

</head>

<body>

<h1>PrimerForge Report</h1>

<p>
Generated:
{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
</p>

<p>
Primer pairs found:
<b>{len(pairs)}</b>
</p>

<table>

<tr>
<th>#</th>
<th>Forward Primer</th>
<th>Reverse Primer</th>
<th>Product (bp)</th>
<th>Forward Tm</th>
<th>Reverse Tm</th>
<th>Forward GC%</th>
<th>Reverse GC%</th>
<th>Population Conservation</th>
<th>Population Coverage</th>
<th>Specificity Score</th>
<th>Specificity</th>
<th>Overall Score</th>
</tr>

{rows}

</table>

</body>

</html>
"""

        output.write_text(
            html,
            encoding="utf-8",
        )

        return output

###############################################################################
# NAME: pyp_reports_templates.py
# VERSION: see PyPedal.__version__
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
# Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle,
# 2025-2026. See CHANGELOG.md for a summary of changes.
# SPDX-License-Identifier: LGPL-2.1-or-later
###############################################################################
# FUNCTIONS:
###############################################################################

_pdfSettings = {
    "_pdfCalcs": {"_page_width": 612, "_page_height": 792},
    "_pdfTitle": "PyPedal PDF Report",
    "_pdfPageinfo": "PyPedal",
}

try:
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import inch
    styles = getSampleStyleSheet()
except ImportError:
    SimpleDocTemplate = Paragraph = Spacer = inch = None
    styles = None

##
# Define the fixed features of the first page of the document with this function.
def myFirstPage(_pdfSettings, canvas, doc):
    """
    Define the fixed features of the first page of the document.
    """
    canvas.saveState()
    canvas.setFont('Times-Bold', 16)
    canvas.drawCentredString(
        _pdfSettings['_pdfCalcs']['_page_width'] / 2.0,
        _pdfSettings['_pdfCalcs']['_page_height'] - 108,
        _pdfSettings['_pdfTitle']
    )
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(inch, 0.75 * inch, f"First Page / {_pdfSettings['_pdfPageinfo']}")
    canvas.restoreState()

##
# Define an alternate layout for the fixed features of subsequent pages.
def myLaterPages(_pdfSettings, canvas, doc):
    """
    Define an alternate layout for the fixed features of subsequent pages.
    """
    canvas.saveState()
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(inch, 0.75 * inch, f"Page {doc.page} {_pdfSettings['_pdfPageinfo']}")
    canvas.restoreState()

##
# Generate the PDF using the settings provided.
def go(_pdfSettings):
    if SimpleDocTemplate is None:
        raise ImportError("ReportLab is required for PDF reports. Install with: pip install 'PyPedal[reports]'")
    output_file = "phello.pdf"
    doc = SimpleDocTemplate(output_file)
    print(f"Writing PDF to {output_file}")
    Story = [Spacer(1, 2 * inch)]
    style = styles["Normal"]

    for i in range(100):
        bogustext = (f"Paragraph number {i}. " * 20)
        p = Paragraph(bogustext, style)
        Story.append(p)
        Story.append(Spacer(1, 0.2 * inch))

    doc.build(
        Story,
        onFirstPage=lambda canvas, doc: myFirstPage(_pdfSettings, canvas, doc),
        onLaterPages=lambda canvas, doc: myLaterPages(_pdfSettings, canvas, doc)
    )

if __name__ == "__main__":
    go(_pdfSettings)
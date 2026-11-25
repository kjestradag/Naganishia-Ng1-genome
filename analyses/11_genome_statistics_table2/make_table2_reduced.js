const fs = require('fs');
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType, AlignmentType,
        BorderStyle, PageOrientation, VerticalAlign } = require('docx');

// [italic, plain, accession, size, max, N50, L50, N90kb, L90, scaff, GC, CDS, BUSCO, bold]
const rows = [
 ["N. adeliensis", " Ng1", "This study (PRJNA1274400)", "19.69","2.25","1.16","7","433.84","18","40","0.53","6,818ᵃ","92.2ᵇ", true],
 ["N. uzbekistanensis", "", "GCA_042257065.1", "24.37","1.79","1.25","8","418.74","22","35","0.53","–","93.9"],
 ["N. vishniacii", " ANT03-052", "GCA_015708705.1", "19.68","1.84","1.08","8","412.92","19","38","0.53","–","93.3"],
 ["N. liquefaciens", "", "GCA_013423385.1", "19.49","1.63","1.03","8","436.83","20","59","0.53","6,883","92.1"],
 ["N. diffluens", "", "GCA_037832395.1", "19.09","1.79","0.80","8","320.67","21","40","0.54","–","92.3"],
 ["N. albida", "", "GCA_001599735.1", "20.70","1.58","0.98","9","241.93","24","74","0.54","–","92.1"],
 ["N. vaughanmartiniae", "", "GCA_030039245.1", "23.00","2.05","0.89","9","255.01","24","64","0.51","6,621","93.1"],
 ["Naganishia", " sp. PG3902-1", "GCA_023628075.1", "19.04","1.38","0.66","9","150.83","27","84","0.54","–","92.1"],
 ["N. friedmannii", "", "GCA_030039255.1", "23.84","2.01","0.75","12","152.95","37","102","0.51","7,205","93.7"],
 ["N. onofrii", "", "GCA_030039215.1", "23.78","2.13","0.69","12","218.90","33","68","0.51","6,962","93.5"],
 ["N. adeliensis", " IF1SW-F1", "GCA_012922715.1", "19.40","1.17","0.51","13","104.16","43","149","0.53","–","92.3", true],
 ["Naganishia", " sp. IF7SW-B1", "GCA_012922685.1", "19.43","1.19","0.51","13","113.41","40","123","0.53","–","92.1"],
 ["N. randhawae", "", "GCA_013461525.1", "20.27","0.79","0.19","30","31.86","134","386","0.52","–","93.1"],
 ["Naganishia", " sp. PT2105A-2", "GCA_023628115.1", "19.39","0.58","0.19","32","55.90","102","228","0.54","–","89.3"],
 ["Naganishia", " sp. PT2105A", "GCA_023628055.1", "18.63","0.53","0.16","36","57.00","110","209","0.54","–","88.7"],
 ["N. kalamii", "", "GCA_022813745.1", "27.86","1.00","0.19","38","30.20","170","478","0.54","–","92.8"],
 ["N. antarctica", "", "GCA_025758675.1", "21.22","0.51","0.17","42","36.00","146","401","0.52","–","92.4"],
 ["N. cerealis", "", "GCA_030039265.1", "24.98","0.56","0.18","49","84.60","130","165","0.49","9,284","83.7"],
];
// reduced version: without N90, scaffold number and GC
const COLS = [3, 4, 5, 6, 8, 11, 12];   // row indices kept: size, largest, N50, L50, L90, CDS, BUSCO
const head = ["Organism","Accession no.","Assembly size (Mb)","Largest scaffold (Mb)","N50 (Mb)","L50","L90","CDS","BUSCO (%)"];
const widths = [3100, 3100, 1250, 1300, 1100, 900, 900, 1150, 1100];
const FONT = "Arial", SZ = 17;
const none = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const line = (sz) => ({ style: BorderStyle.SINGLE, size: sz, color: "000000" });

function cell(children, w, opts = {}) {
  return new TableCell({
    width: { size: w, type: WidthType.DXA }, verticalAlign: opts.va || VerticalAlign.CENTER,
    margins: { top: 30, bottom: 30, left: 70, right: 70 },
    borders: { top: opts.top || none, bottom: opts.bottom || none, left: none, right: none },
    children: [new Paragraph({ alignment: opts.align || AlignmentType.LEFT, spacing: { before: 0, after: 0 }, children })],
  });
}
const hdr = new TableRow({ tableHeader: true, children: head.map((h, i) =>
  cell([new TextRun({ text: h, bold: true, font: FONT, size: SZ })], widths[i],
       { top: line(12), bottom: line(6), align: i < 2 ? AlignmentType.LEFT : AlignmentType.RIGHT, va: VerticalAlign.BOTTOM })) });
const body = rows.map((r, k) => {
  const last = k === rows.length - 1, b = !!r[13];
  const brd = last ? { bottom: line(12) } : {};
  const org = [new TextRun({ text: r[0], italics: true, bold: b, font: FONT, size: SZ })];
  if (r[1]) org.push(new TextRun({ text: r[1], bold: b, font: FONT, size: SZ }));
  const cells = [cell(org, widths[0], brd), cell([new TextRun({ text: r[2], bold: b, font: FONT, size: SZ })], widths[1], brd)];
  COLS.forEach((i, j) =>
    cells.push(cell([new TextRun({ text: r[i], bold: b, font: FONT, size: SZ })], widths[j + 2], { ...brd, align: AlignmentType.RIGHT })));
  return new TableRow({ children: cells });
});
const total = widths.reduce((a, b) => a + b, 0);
const table = new Table({ width: { size: total, type: WidthType.DXA }, columnWidths: widths, rows: [hdr, ...body] });
const P = (runs, after = 60) => new Paragraph({ spacing: { before: 0, after }, children: runs });
const t = (text, o = {}) => new TextRun({ text, font: FONT, size: o.size || 18, bold: o.bold, italics: o.italics });
const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: 18 } } } },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838, orientation: PageOrientation.LANDSCAPE },
                          margin: { top: 1134, bottom: 1134, left: 1134, right: 1134 } } },
    children: [
      P([t("Table 2. ", { bold: true }), t("Assembly statistics of "), t("N. adeliensis", { italics: true }),
         t(" Ng1 and publicly available "), t("Naganishia", { italics: true }), t(" genomes, ordered by L50 (ties broken by N50).")], 120),
      table,
      P([], 60),
      P([t("BUSCO, complete single-copy and duplicated genes (BUSCO v5.8.0, basidiomycota_odb10, genome mode with miniprot, same settings for all assemblies); CDS, protein-coding genes in the NCBI annotation; –, no annotation available. Rows in bold: strains of ", { size: 16 }),
         t("N. adeliensis", { italics: true, size: 16 }), t(".", { size: 16 })]),
      P([t("ᵃ 6,774 nuclear and 44 mitochondrial protein-coding genes.", { size: 16 })]),
      P([t("ᵇ 93.4% on the annotated Ng1 proteome (protein mode).", { size: 16 })]),
    ],
  }],
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(process.argv[2], b); console.log("written", process.argv[2], "table width", total); });

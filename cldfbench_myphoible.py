"""
Generates a CLDF dataset with data derived from Moran & McCloy's "Phoible 2.0" (2020).

This complex method is only needed while we cannot get raw strings (pre-processing,
pre-normalization) from Phoible.
"""

# TODO: do allophones, which are currently just stripped off

import csv
import unicodedata
from pathlib import Path
from unidecode import unidecode
from collections import defaultdict

from pyglottolog import Glottolog
from pyclts import CLTS, models

from pycldf import Sources
from cldfbench import CLDFSpec
from cldfbench import Dataset as BaseDataset
from clldutils.misc import slug


def compute_id(text):
    """
    Returns a codepoint representation to an Unicode string.
    """

    unicode_repr = "".join(["u{0:0{1}X}".format(ord(char), 4) for char in text])

    label = slug(unidecode(text))

    return "%s_%s" % (label, unicode_repr)


def normalize_grapheme(text):
    """
    Apply simple, non-CLTS, normalization.
    """

    text = unicodedata.normalize("NFC", text)

    if text[0] == "(" and text[-1] == ")":
        text = text[1:-1]

    if text[0] == "[" and text[-1] == "]":
        text = text[1:-1]

    return text


def read_phoible_common(filename, source):
    """
    Reads common Phoible sources.

    Sources that can be read with this function include AA, EA, ER, GM, PH, and UZ.
    """

    # Read lines one by one, fixing issues in the TSV file
    curr_inventory_id = None
    curr_lang_name = None
    segment_field = False

    entries = []
    with open(filename) as tsvfile:
        # Build reader and set the raw source column
        reader = csv.DictReader(tsvfile, delimiter="\t")
        if "PhonemeOld" in reader.fieldnames:
            segment_field = "PhonemeOld"
        elif "Segment" in reader.fieldnames:
            segment_field = "Segment"
        else:
            segment_field = "Phoneme"

        for row in reader:
            # Store the current InventoryID and LanguageName, if we have a new one (it
            # is only given for the first row of each inventory)
            if row["InventoryID"]:
                curr_inventory_id = row["InventoryID"]
            if row["LanguageName"]:
                curr_lang_name = row["LanguageName"]

            # Build a new entry
            segment_raw = row[segment_field]
            if segment_raw:
                entries.append(
                    {
                        "lang_id": f"{source}_{curr_lang_name}",
                        "source": f"PHOIBLE_{source}",
                        "inv_id": curr_inventory_id,
                        "segment_raw": segment_raw,
                    }
                )

    return entries


def read_phoible_ra(filename):
    """
    Reads common Phoible RA sources.
    """

    entries = []
    with open(filename) as csvfile:
        # The first row has full segments description, the second and third have the
        # actual headers
        reader = csv.reader(csvfile)
        next(reader)
        next(reader)

        # Get field name and full segment list
        fields = next(reader)

        # Iterate over data
        for row in reader:
            data = {field: value for field, value in zip(fields, row)}

            inventory_id = data["InventoryID"]
            language_name = data["Language Name"]

            for segment_raw in fields:
                if segment_raw not in [
                    "InventoryID",
                    "#",
                    "Language Name",
                    "Language Code",
                ]:
                    if data[segment_raw] == "1":
                        entry = {
                            "lang_id": f"RA_{language_name}",
                            "source": "PHOIBLE_RA",
                            "inv_id": inventory_id,
                            "segment_raw": segment_raw,
                        }

                        # Update entries
                        entries.append(entry)

    return entries


def read_phoible_saphon(filename):
    """
    Reads common Phoible SAPHON sources.
    """

    entries = []
    with open(filename) as tsvfile:
        # Build reader and get list of segments
        reader = csv.DictReader(tsvfile, delimiter="\t")
        segments = reader.fieldnames[17:-37]

        for row in reader:
            for segment_raw in segments:
                if row[segment_raw] == "1":

                    entry = {
                        "lang_id": f"SAPHON_{row['Name']}",
                        "source": "PHOIBLE_SAPHON",
                        "inv_id": row["InventoryID"],
                        "segment_raw": segment_raw,
                    }

                    # Update entries
                    entries.append(entry)

    return entries


def read_phoible_upsid(dirname):
    """
    Reads common Phoible UPSID sources.
    """

    entries = []

    # load list of languages
    lang_name2num = {}
    lang_num2name = {}
    with open(dirname / "UPSID_Languages.tsv") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            lang_name2num[row["LangName"]] = {"num": row["LangNum"]}
            lang_num2name[row["LangNum"]] = {"name": row["LangName"]}

    # load inventory id
    with open(dirname / "UPSID_LanguageCodes.tsv") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            lang_name2num[row["LanguageName"]]["id"] = row["InventoryID"]
            lang_num2name[row["upsidLangNum"]]["id"] = row["InventoryID"]

    # load UPSID CC ID
    # TODO: study why Python csv is failing
    segments = {}
    with open(dirname / "UPSID_CharCodes.tsv") as tsvfile:
        fields = next(tsvfile).split("\t")
        for row in tsvfile:
            values = row.split("\t")
            data = {field: value for field, value in zip(fields, values)}

            segments[data["CCID"]] = data["IPA"]

    # load segments
    with open(dirname / "UPSID_Segments.tsv") as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        for row in reader:
            inv_id = lang_num2name[row["upsidLangNum"]]["id"]
            entry = {
                "lang_id": f"UPSID_{lang_num2name[row['upsidLangNum']]['name']}",
                "source": "PHOIBLE_UPSID",
                "inv_id": inv_id,
                "segment_raw": segments[row["upsidCCID"]],
                "marginal": str(row["anomalous"] != "0"),
            }

            entries.append(entry)

    return entries


def read_raw_phoible_data(raw_dir):
    data = []

    # Read Phoible AA data
    filename = raw_dir / "phoible-dev" / "raw-data" / "AA" / "AA_inventories.tsv"
    data += read_phoible_common(filename, "AA")

    # Read Phoible EA data
    filename = raw_dir / "phoible-dev" / "raw-data" / "EA" / "EA_inventories.tsv"
    data += read_phoible_common(filename, "EA")

    # Read Phoible ER data
    filename = raw_dir / "phoible-dev" / "raw-data" / "ER" / "ER_inventories.tsv"
    data += read_phoible_common(filename, "ER")

    # Read Phoible GR-afr data
    filename = raw_dir / "phoible-dev" / "raw-data" / "GM" / "gm-afr-inventories.tsv"
    data += read_phoible_common(filename, "GM")

    # Read Phoible GR-sea data
    filename = raw_dir / "phoible-dev" / "raw-data" / "GM" / "gm-sea-inventories.tsv"
    data += read_phoible_common(filename, "GM")

    # Read Phoible PH data
    filename = raw_dir / "phoible-dev" / "raw-data" / "PH" / "phoible_inventories.tsv"
    data += read_phoible_common(filename, "PH")

    # Read Phoible UW data
    filename = raw_dir / "phoible-dev" / "raw-data" / "UZ" / "UZ_inventories.tsv"
    data += read_phoible_common(filename, "UZ")

    # Read Phoible RA data
    filename = raw_dir / "phoible-dev" / "raw-data" / "RA" / "Ramaswami1999.csv"
    data += read_phoible_ra(filename)

    # Read Phoible SAPHON data
    filename = raw_dir / "phoible-dev" / "raw-data" / "SAPHON" / "saphon20121031.tsv"
    data += read_phoible_saphon(filename)

    # Read Phoible UPSID data
    dirname = raw_dir / "phoible-dev" / "raw-data" / "UPSID"
    data += read_phoible_upsid(dirname)

    return data


class Dataset(BaseDataset):
    """
    CLDF dataset for inventories.
    """

    dir = Path(__file__).parent
    id = "myphoible"

    def cldf_specs(self):  # A dataset must declare all CLDF sets it creates.
        return CLDFSpec(dir=self.cldf_dir, module="StructureDataset")

    def cmd_download(self, args):
        """
        Download files to the raw/ directory. You can use helpers methods of `self.raw_dir`, e.g.

        >>> self.raw_dir.download(url, fname)
        """
        pass

    def cmd_makecldf(self, args):
        """
        Convert the raw data to a CLDF dataset.

        >>> args.writer.objects['LanguageTable'].append(...)
        """

        # Instantiate Glottolog and CLTS
        # TODO: how to call CLTS?
        glottolog = Glottolog(args.glottolog.dir)
        clts_path = Path.home() / ".config" / "cldf" / "clts"
        clts = CLTS(clts_path)

        # Add the bibliographic info, loading sources and bibliographic info
        # provided in phoible/dev
        source_map = defaultdict(list)
        with open(
            self.raw_dir / "phoible-dev" / "mappings" / "InventoryID-Bibtex.csv"
        ) as csvfile:
            for row in csv.DictReader(csvfile):
                source_map[row["InventoryID"]].append(row["BibtexKey"])
        sources = Sources.from_file(
            self.raw_dir / "phoible-dev" / "data" / "phoible-references.bib"
        )
        args.writer.cldf.add_sources(*sources)

        # Add components
        args.writer.cldf.add_columns(
            "ValueTable",
            {"name": "Marginal", "datatype": "boolean"},
            "Catalog",
            "Contribution_ID",
        )

        args.writer.cldf.add_component("ParameterTable", "BIPA")
        args.writer.cldf.add_component(
            "LanguageTable", "Family_Glottocode", "Family_Name", "Glottolog_Name"
        )
        args.writer.cldf.add_table(
            "inventories.csv",
            "ID",
            "Name",
            "Contributor_ID",
            {
                "name": "Source",
                "propertyUrl": "http://cldf.clld.org/v1.0/terms.rdf#source",
                "separator": ";",
            },
            "URL",
            "Tones",
            primaryKey="ID",
        )

        # Read Phoible language mapping
        phoible_lang = {}
        glottocodes = set()
        for row in self.etc_dir.read_csv("languages.csv", dicts=True):
            # inventory to language mapping; remember that in many cases we have
            # more than one inventory per glottocode
            inv_id = row["ID"].split("_")[0]
            phoible_lang[inv_id] = row["ID"]

            # We collect all glottocodes, so later we can quicker access the Glottolog
            # data without traversing N times
            if row["Glottocode"]:
                glottocodes.add(row["Glottocode"])

        # Build language info, after caching glottolog langouids
        languages = []
        languoids = {}

        for languoid in glottolog.languoids():
            if languoid.glottocode in glottocodes:
                languoids[languoid.glottocode] = languoid

        for row in self.etc_dir.read_csv("languages.csv", dicts=True):
            lang_dict = {"ID": row["ID"], "Name": row["Name"]}
            if row["Glottocode"]:
                lang = languoids[row["Glottocode"]]
                lang_dict.update(
                    {
                        "Family_Glottocode": lang.lineage[0][1]
                        if lang.lineage
                        else None,
                        "Family_Name": lang.lineage[0][0] if lang.lineage else None,
                        "Glottocode": row["Glottocode"],
                        "ISO639P3code": lang.iso_code,
                        "Latitude": lang.latitude,
                        "Longitude": lang.longitude,
                        "Macroarea": lang.macroareas[0].name
                        if lang.macroareas
                        else None,
                        "Glottolog_Name": lang.name,
                    }
                )
            languages.append(lang_dict)

        # Read Phoible data
        data = read_raw_phoible_data(self.raw_dir.absolute())

        # Add data to CLDF
        parameters = []
        values = []
        counter = 1
        for entry in data:
            # Mark marginality, if necessary
            if entry["segment_raw"][0] == "<" and entry["segment_raw"][-1] == ">":
                entry["segment_raw"] = entry["segment_raw"][1:-1]
                entry["marginal"] = "True"
            elif entry["segment_raw"][0] == "'" and entry["segment_raw"][-1] == "'":
                entry["segment_raw"] = entry["segment_raw"][1:-1]
                entry["marginal"] = "True"
            else:
                entry["marginal"] = "False"

            # Use the first phoneme of an allophone group, if any
            if len(entry["segment_raw"]) >= 3 and "|" in entry["segment_raw"]:
                entry["segment_raw"] = entry["segment_raw"].split("|")[0]

            # Obtain the corresponding BIPA grapheme, is possible
            normalized = normalize_grapheme(entry["segment_raw"])
            sound = clts.bipa[normalized]
            if isinstance(sound, models.UnknownSound):
                par_id = "UNK_" + compute_id(normalized)
                bipa_grapheme = ""
                desc = ""
            else:
                par_id = "BIPA_" + compute_id(normalized)
                bipa_grapheme = str(sound)
                desc = sound.name
            parameters.append((par_id, normalized, bipa_grapheme, desc))

            values.append(
                {
                    "ID": str(counter),
                    "Language_ID": phoible_lang[entry["inv_id"]],
                    "Parameter_ID": par_id,
                    "Value": entry["segment_raw"],
                    "Contribution_ID": entry["inv_id"],
                    "Source": source_map[entry["inv_id"]],
                    "Catalog": slug(entry["source"]),
                }
            )
            counter += 1

        # Build segment data
        segments = [
            {"ID": id, "Name": normalized, "BIPA": bipa_grapheme, "Description": desc}
            for id, normalized, bipa_grapheme, desc in set(parameters)
        ]

        # Write data and validate
        args.writer.write(
            **{
                "ValueTable": values,
                "LanguageTable": languages,
                "ParameterTable": segments,
                "inventories.csv": [],
            }
        )

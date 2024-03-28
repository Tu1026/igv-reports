import os
import sys
import json
import math
import argparse
from urllib.request import urlopen
from igv_reports import datauri, tracks, utils
from igv_reports.varianttable import VariantTable
from igv_reports.bedtable import BedTable
from igv_reports.bedtable import BedpeTable
from igv_reports.bedtable import JunctionBedTable
from igv_reports.generictable import GenericTable
from igv_reports.regions import parse_region
from igv_reports.feature import MockReader
from igv_reports.fasta import FastaReader
from igv_reports.ideogram import IdeogramReader
from igv_reports.genome import get_genome
from igv_reports.tracks import get_track_type, is_format_supported
# import asyncio
# import threading
import sys
import multiprocessing as mp
from itertools import repeat
from pebble import ProcessPool,ProcessExpired
from pathlib import Path
import traceback
### Speed things up with Async
from concurrent.futures import TimeoutError


'''
Create an html report.  This is the main function for the application.
'''
def create_report(args):

    trackconfigs = []
    trackreaders = []
    session_dict = {}

    # Read the variant data
    variants_file = args.sites

    if variants_file.endswith(".bcf") or variants_file.endswith(".vcf") or variants_file.endswith(".vcf.gz"):
        table = VariantTable(variants_file, args.info_columns, args.info_columns_prefixes, args.samples,
                             args.sample_columns, args.idlink)

    elif variants_file.endswith(".bed") or variants_file.endswith(".bed.gz"):
        if args.type is not None and args.type == "junction":
            table = JunctionBedTable(variants_file, args.info_columns)
        else:
            table = BedTable(variants_file)

    elif variants_file.endswith(".bedpe") or variants_file.endswith(".bedpe.gz"):
        table = BedpeTable(variants_file)

    elif variants_file.endswith(".maf") or variants_file.endswith(".maf.gz") or (
            args.sequence is not None and args.begin is not None and args.end is not None):
        # A hack -- since these file formats are not supported by igv.js mock up a bed style track from a table of generic features
        table = GenericTable(variants_file, args.info_columns, args.sequence, args.begin, args.end, args.zero_based, args.bam)
        flist = []
        for tuple in table.features:
            flist.append(tuple[0])
        trackconfigs.append({
            "config": {"name": "variants", "type": "annotation", "format": "bed"},
            "reader": MockReader(flist)
        })

    # Track json array.  Tracks can come from (1) tracks CL argument, (2) genome CL argument, and (3) track_config Cl argument
    trackjson = []

    # Check for optional genome argument.  If supplied an igv.js genome json definition is used in lieu of a fasta file
    if args.genome is not None:
        genome = get_genome(args.genome)
        if args.fasta is None:
            args.fasta = genome["fastaURL"]
        if args.ideogram is None and "cytobandURL" in genome:
            args.ideogram = genome["cytobandURL"]
        if "tracks" in genome:
            for config in genome["tracks"]:
                if "format" not in config and "url in c":
                    config["format"] = feature.infer_format(config["url"])
                if is_format_supported(config["format"]):
                    if "type" not in config:
                        config["type"] = get_track_type(config["format"])
                    trackjson.append(config)
                else:
                    print("File format: " + config["format"] + " is not supported. Skipping track '" + config["name"] + "'.")
                    
    if args.tracks is not None:
        for track in args.tracks:
           trackjson.append(tracks.get_track_json_dict(track))

    if args.track_config is not None:
        for trackobj in args.track_config:
            with open(trackobj) as f:
                j = json.load(f)
                for c in j:
                    trackjson.append(c)


    # Create file readers for tracks.  This is done outside the locus loop so initialization happens once

    for config in trackjson:
        reader = utils.getreader(config, None, args)
        trackconfigs.append({
            "config": config,
            "reader": reader
        })



    # Other readers
    # fasta_reader = FastaReader(args.fasta)
    # if (args.ideogram):
    #     ideogram_reader = IdeogramReader(args.ideogram)
    # else:
    #     ideogram_reader = None

    # loop through regions defined from variant, annotation, or bedpe files,  creating an igv.js session for each one
    flanking = 0
    if args.flanking is not None:
        flanking = float(args.flanking)

    #TODO fix the ugly multiprocessing code
    execute_status= []
    cores = args.n_workers
    TIMEOUT_SECONDS = 60
    n_tasks = len(table.features)
    results = []
    with mp.Manager() as manager:
        sh_session_dict = manager.dict()
        sh_execute_status = manager.list()
        # with manager.Pool(processes=cores) as pool:
        futures = []
        with ProcessPool(max_workers=cores) as pool:
            i = 0
            for tuple in table.features:
                i += 1
                future = pool.schedule(feature_process, args=(i, tuple, sh_session_dict, sh_execute_status, trackconfigs, n_tasks, flanking, args), timeout=60)
                # future.add_done_callback(task_done)
                futures.append(future)
            # future = pool.map(feature_process, [list(range(1, n_tasks+1)), list(table.features), repeat(sh_session_dict), repeat(sh_execute_status), repeat(trackconfigs), repeat(n_tasks), repeat(flanking), repeat(args)], timeout=60)
            # future = pool.startmap(feature_process, zip(list(range(1, n_tasks+1)), list(table.features), repeat(sh_session_dict), repeat(sh_execute_status), repeat(trackconfigs), repeat(n_tasks), repeat(flanking), repeat(args)), timeout=60)

            for i, f in enumerate(futures):
                try:
                    results.append(future.result(TIMEOUT_SECONDS))
                except TimeoutError:
                    execute_status.append([f"{i}th variant failed because it timed out after {TIMEOUT_SECONDS} seconds"])
                except Exception as error:
                    print("Exception type: ", type(error))  # The type of the exception
                    print("Exception args: ", error.args)  # The arguments that the exception has been called with
                    print("Exception message: ", str(error))  # The error message
                    print("Traceback: ")
                    traceback.print_exc()  # Python's traceback of remote process
                    execute_status.append([f"{i}th variant failed because of an error", str(table.features[i][0])])

                # execute_status.extend(list(sh_execute_status))
                # session_dict = {**session_dict, **dict(sh_session_dict)}
                
        execute_status.extend(list(sh_execute_status))
        session_dict = {**session_dict, **dict(sh_session_dict)}
        print(results)
    # with mp.Manager() as manager:
    #     sh_session_dict = manager.dict()
    #     sh_execute_status = manager.list()
    #     with manager.Pool(processes=cores) as pool:
    #     # with ProcessPool(max_workers=cores) as pool:
    #         i = 0
    #         futures = []
    #         results = []
    #         timed_out_results=[]
    #         for tuple in table.features:
    #             i += 1
    #             # futures.append(pool.schedule(feature_process, args=(i, tuple, sh_session_dict, sh_execute_status, trackconfigs, n_tasks, flanking, args), timeout=60))
    #             futures.append(pool.apply_async(feature_process, args=(i, tuple, sh_session_dict, sh_execute_status, trackconfigs, n_tasks, flanking, args)))
    #         # future = pool.map(feature_process, zip(list(range(1, n_tasks+1)), list(table.features), repeat(sh_session_dict), repeat(sh_execute_status), repeat(trackconfigs), repeat(n_tasks), repeat(flanking), repeat(args)), timeout=60)
    #         # future = pool.startmap(feature_process, zip(list(range(1, n_tasks+1)), list(table.features), repeat(sh_session_dict), repeat(sh_execute_status), repeat(trackconfigs), repeat(n_tasks), repeat(flanking), repeat(args)), timeout=60)
    #         for i, future in enumerate(futures):
    #             try:
    #                 results.append(future.get(60))
    #             except TimeoutError:
    #                 timed_out_results.append([future, f"{i}th variant failed"])

    #             execute_status.extend(list(sh_execute_status))
    #             session_dict = {**session_dict, **dict(sh_session_dict)}
    #     print(results)
        
    output_folder = Path(args.output.split(".")[0])
    
    ## Creat folder.    
    output_folder.mkdir(parents=False, exist_ok=True)
    
    ## Add error tracking functionality
    with open(output_folder.joinpath("error_files.txt"), "w") as f:
        for status in execute_status:
            if status[0] == 1:
                # f.write(f"Error: {status[1]} \n")
                f.write(f"Traceback: {status[2]} \n")
                f.write(f"Variant: {status[3]} \n")
                f.write(f"Failed_feature: {str(status[4])} \n")
                f.write("\n")
            
    # session_dict = json.dumps(session_dict)

    
    print("start dumping")
    key_dict, first_chunk, first_file = chunk_dict(session_dict, 200, output_folder)
    session_dict = first_chunk
    print("finished dumping")

    template_file = args.template
    if None == template_file:
        if 'junction' == args.type:
            template_file = os.path.dirname(sys.modules['igv_reports'].__file__) + '/templates/junction_template.html'
        else:
            template_file = os.path.dirname(sys.modules['igv_reports'].__file__) + '/templates/variant_template-2.html'

    output_file = args.output

    standalone = args.standalone
    with open(template_file, "r") as f:
        data = f.readlines()

        with open(output_folder.joinpath(output_file), "w") as o:

            for i, line in enumerate(data):

                if args.title is not None and line.startswith("<!--title-->"):
                    o.write("<h1>" + args.title + "</h1>")

                if standalone:
                    if line.strip().startswith("<script") and ".js\"" in line:
                        inline_script(line, o, "js")
                        continue
                    elif line.strip().startswith("<link") and line.strip().endswith("css\">"):
                        inline_script(line, o, "css")
                        continue
                j = line.find('"@TABLE_JSON@"')
                if j >= 0:
                    line = line.replace('"@TABLE_JSON@"', table.to_JSON())

                j = line.find('"@SESSION_DICTIONARY@"')
                if j >= 0:
                    line = line.replace('"@SESSION_DICTIONARY@"', session_dict)

                j = line.find('"@KEY_DICTIONARY@"')
                if j >= 0:
                    line = line.replace('"@KEY_DICTIONARY@"', key_dict)  
                j = line.find('"@CURR_FILE@"')
                if j >= 0:
                    line = line.replace('"@CURR_FILE@"', first_file)  
                o.write(line)

def chunk_dict(session_dict, chunk_size, folder):
    chunks = [dict(list(session_dict.items())[i:i + chunk_size]) for i in range(0, len(session_dict), chunk_size)]
    key_to_file = {}

    for i, chunk in enumerate(chunks):
        filename = f'chunk{i}.json'
        with open(folder.joinpath(filename), 'w') as f:
            json.dump(chunk, f)
        for key in chunk.keys():
            key_to_file[key] = filename

        
    return  json.dumps(key_to_file), json.dumps(chunks[0]), f'"chunk{0}.json"'


def task_done(future):
    try:
        result = future.result()  # blocks until results are ready
    except TimeoutError as error:
        print("Function took longer than %d seconds" % error.args[1])
    except Exception as error:
        print(f"Function raised {error}")
        error.__traceback__.print_tb()  # traceback of the function
        
def inline_script(line, o, source_type):
    # <script type="text/javascript" src="https://igv.org/web/test/dist/igv.min.js"></script>
    if source_type == "js":
        s = line.find('src="')
        offset = 5
        o.write('<script type="text/javascript">\n')
    elif source_type == "css":
        s = line.find('href="')
        offset = 6
        o.write('<style type="text/css">\n')
    else:
        raise KeyError("Inline script must be either js- or css-file")
    if s > 0:
        e = line.find('">', s)
        url = line[s + offset:e]
        response = urlopen(url)
        content = response.read().decode('utf-8')
        response.close()
        o.write(content)
        if source_type == "js":
            o.write('</script>\n')
        else:
            o.write('</style>\n')
    else:
        raise ValueError("No file path in {l} for inline script.".format(l=line))


def locus_string(chr, start, end):
    if (end - start) == 1:
        return f'{chr}:{start + 1}'
    else:
        return f'{chr}:{start + 1}-{end}'


def feature_process(i, tuple,session_dict, sh_execute_status,trackconfigs, n_tasks, flanking, args):
    # Other readers
    fasta_reader = FastaReader(args.fasta)
    if (args.ideogram):
        ideogram_reader = IdeogramReader(args.ideogram)
    else:
        ideogram_reader = None
    result = [0, ""]
    try: 
        print(f"Working on variant {i}/{n_tasks}")

        feature = tuple[0]
        unique_id = tuple[1]

        # If a variant feature (=> row in the table) has an explicit session id use it, otherwise use the row id (i.e. unique_id)
        if hasattr(feature, "session_id"):
            session_id = feature.session_id
        else:
            session_id = str(unique_id)
            

        if session_id not in session_dict:

            # Placeholder variable for possible second region (bedpe files)
            region2 = None

            # Define a genomic region around the variant
            if hasattr(feature, "viewport"):
                region = parse_region(feature.viewport)
                chr = region["chr"]
                start = region["start"]
                end = region["end"]
            else:
                chr = feature.chr
                start = int(math.floor(feature.start - flanking / 2))
                start = max(start, 1)  # bound start to 1
                end = int(math.ceil(feature.end + flanking / 2))
                region = {
                    "chr": chr,
                    "start": start,
                    "end": end
                }

                # If feature has a second locus (bedpe file) create the region here.
                if hasattr(feature, 'chr2') and feature.chr2 is not None:
                    chr2 = feature.chr2
                    start2 = int(math.floor(feature.start2 - flanking / 2))
                    start2 = max(start2, 1)  # bound start to 1
                    end2 = int(math.ceil(feature.end2 + flanking / 2))
                    region2 = {
                        "chr": chr2,
                        "start": start2,
                        "end": end2
                    }

            data = fasta_reader.slice(region)
            fa = '>' + chr + ':' + str(start) + '-' + str(end) + '\n' + data

            if region2 is not None:
                data2 = fasta_reader.slice(region2)
                fa += '\n' + '>' + chr2 + ':' + str(start2) + '-' + str(end2) + '\n' + data2

            fasta_uri = datauri.get_data_uri(fa)
            fastaJson = {
                "fastaURL": fasta_uri,
            }

            # Ideogram
            if (args.ideogram):
                ideo_string = ideogram_reader.get_data(region["chr"])
                if region2 is not None:
                    ideo_string += ideogram_reader.get_data(region2["chr"])
                ideo_uri = datauri.get_data_uri(ideo_string)
                fastaJson["cytobandURL"] = ideo_uri

            # Initial locus
            if (hasattr(feature, "viewport")):
                initial_locus = feature.viewport
            else:
                initial_locus = locus_string(feature.chr, feature.start, feature.end)
                if region2 is not None:
                    initial_locus += f' {locus_string(feature.chr2, feature.start2, feature.end2)}'

            session_json = {
                "locus": initial_locus,
                "reference": fastaJson,
                "tracks": []
            }
            
            trackconf = []
            if args.bam is not None:
                config = tracks.get_track_json_dict(feature.bam)
                reader = utils.getreader(config, None, args)
                trackconf.append({
                    "config": config,
                    "reader": reader
                })

            track_objects = []
            # Loop through user supplied track configs
            # "cram" input format is converted to "bam" for output track configs
            for tc in [*trackconfigs, *trackconf]:
                trackobj = tc["config"];

                # Set some defaults if not specified
                if "url" in trackobj:
                    default_trackobj = tracks.get_track_json_dict(trackobj["url"]);
                    if "type" not in trackobj:
                        trackobj["type"] = default_trackobj["type"]
                    if "format" not in trackobj:
                        trackobj["format"] = default_trackobj["format"]
                    if trackobj["format"] == "cram":
                        trackobj["format"] = "bam"
                    if "name" not in trackobj:
                        trackobj["name"] = default_trackobj["url"]

                # Indexes are not used with data URIs
                if "indexURL" in trackobj:
                    del trackobj["indexURL"]

                reader = tc["reader"]
                data = reader.slice(region, region2)
                trackobj["url"] = datauri.get_data_uri(data)
                track_objects.append(trackobj)

            track_order = 1
            for trackobj in track_objects:
                if (trackobj["type"] == "alignment"):
                    trackobj["height"] = 500
                    is_snv = feature.end - feature.start == 1
                    if (trackobj["type"]) == "alignment" and (args.sort is not None or is_snv) and (
                            args.sort != 'NONE'):
                        sort_option = 'BASE' if args.sort is None else args.sort.upper()
                        trackobj["sort"] = {
                            "option": sort_option,
                            "chr": chr,
                            "position": str(feature.start + 1),
                            "direction": "ASC"
                        }
                if "order" not in trackobj:
                    trackobj["order"] = track_order
                session_json["tracks"].append(trackobj)
                track_order += 1

            # Build the session data URI
            session_string = json.dumps(session_json)
            session_uri = datauri.get_data_uri(session_string)
            session_dict[session_id] = session_uri

            
            
    except Exception as e:
        print("Error in feature_process")
        import traceback
        result = [1, str(e), traceback.format_exc(), f"Variant {i}/{n_tasks} failed", feature]
    
    print(f"Done with variant {i}/{n_tasks}")

    sh_execute_status.append(result)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sites", help="vcf file defining variants, required")
    parser.add_argument("fasta", nargs="?", default=None, help="reference fasta file, required if --genome is not specified")
    parser.add_argument("--genome",  help="igv.js genome id (e.g. hg38)")
    parser.add_argument("--type", help="Report type.  Possible values are mutation and junctions.  Default is mutation")
    parser.add_argument("--ideogram", help="ideogram file in UCSC cytoIdeo format")
    parser.add_argument("--tracks", nargs="+", help="list of track files")
    parser.add_argument("--track-config", nargs="+", help="track json file")
    parser.add_argument("--sort",
                        help="initial sort option for alignment tracks.   Supported values include  BASE, STRAND, INSERT_SIZE, and MATE_CHR. Default value is BASE for single nucleotide variants, no sorting otherwise.  See the igv.js documentation for more information. ")
    parser.add_argument("--template", help="html template file", default=None)
    parser.add_argument("--output", help="output file name", default="igvjs_viewer.html")
    parser.add_argument("--info-columns", nargs="+", help="list of VCF info field names to include in variant table")
    parser.add_argument("--info-columns-prefixes", nargs="+",
                        help="list of prefixes of VCF info field names to include in variant table")
    parser.add_argument("--samples", nargs="+",
                        help="Space delimited list of sample (i.e. genotypes) names.  Used in conjunction with --sample-columns")
    parser.add_argument("--sample-columns", nargs="+",
                        help="list of VCF sample (genomtype) FORMAT field names to include in variant table")
    parser.add_argument("--flanking", help="genomic region to include either side of variant", default=1000)
    parser.add_argument("--standalone", help="embed javascript as well as data in output html", action='store_true')
    parser.add_argument("--title", help="optional title string")
    parser.add_argument("--sequence", help="Column of sequence (chromosome) name.  For tab-delimited sites file.",
                        default=None)
    parser.add_argument("--begin", help="Column of start position.  For tab-delimited sites file.", default=None)
    parser.add_argument("--end", help="Column of end position. For tab-delimited sites file.", default=None)
    parser.add_argument("--bam", type=str, help="Column of bam file uri. For tab-delimited sites file.", default=None)
    parser.add_argument("--zero_based",
                        help="Specify that the position in the data file is 0-based (e.g. UCSC files) rather than 1-based.",
                        default=None)
    parser.add_argument("--idlink", type=str, help="url link template for the VCF ID column")
    parser.add_argument("--exclude-flags", type=int, help="Passed to samtools to filter alignments.  For BAM and CRAM files.", default=1536)
    parser.add_argument("--n_workers", type=int, help="Numer of processes to work on at the same time (use something like 100-200 on the server if you are in a rush)", default=12)
    args = parser.parse_args()
    create_report(args)


if __name__ == "__main__":
    main()

import contextlib
import csv
import datetime
import os
import subprocess
import tarfile
import urllib.error
import urllib.parse
import urllib.request
import shutil
import glob

from retrying import retry

csv.field_size_limit(10000000)

#######################
# Convenience methods #
#######################


def create_path_if_nonexistent(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def run_process(args, redirect_stdout_path=None, expected_returncode=0):
    if redirect_stdout_path:
        stdout_cm = open(redirect_stdout_path, 'w')
    else:
        # wrapping into a dummy context, see https://docs.python.org/3/library/contextlib.html#contextlib.nullcontext
        stdout_cm = contextlib.nullcontext(subprocess.PIPE)

    with stdout_cm as stdout:
        print(f"Running with the following args: {args}")
        sp = subprocess.Popen(args, stdout=stdout, stderr=subprocess.PIPE)
        print_subprocess_output_and_error(sp)

        if sp.returncode != expected_returncode:
            raise RuntimeError(f"Process returned with code {sp.returncode} instead of {expected_returncode}")


def print_subprocess_output_and_error(sp):
    out, err = sp.communicate()
    if out:
        print("standard output of subprocess:")
        print(out)
    if err:
        print("standard error of subprocess:")
        print(err)


@retry(stop_max_attempt_number=3, wait_fixed=3000)
def urlopen_with_retry(url):
    return urllib.request.urlopen(url)


def download_file_and_display_progress(url, file_name=None):
    if file_name is None:
        file_name = url.split('/')[-1]

    u = urlopen_with_retry(url)
    f = open(file_name, 'wb')
    meta = u.info()
    file_size = int(meta.get("Content-Length")[0])
    print("Downloading: %s Bytes: %s" % (file_name, file_size))

    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break

        file_size_dl += len(buffer)
        f.write(buffer)
        status = r"%10d  [%3.2f%%]" % (
        file_size_dl, file_size_dl * 100. / file_size)
        status = status + chr(8) * (len(status) + 1)
        print(status, end=' ')

    f.close()
    print("Finished downloading %s" % (file_name))


def download_file_with_basic_auth(url, file_name, username, password):
    p = urllib.request.HTTPPasswordMgrWithDefaultRealm()

    p.add_password(None, url, username, password)

    handler = urllib.request.HTTPBasicAuthHandler(p)
    opener = urllib.request.build_opener(handler)
    urllib.request.install_opener(opener)

    data = urlopen_with_retry(url).read()
    f = open(file_name, "wb")
    f.write(data)
    f.close()
    print("Finished downloading %s" % (file_name))


def check_file_for_contents(file_path):
    handle_process_success_or_failure(os.stat(file_path).st_size != 0,
                                      file_path)


def concatenate_symbols(symbols):
    concatted_symbols = []
    for symbol in symbols:
        concatted_symbols.append(symbol)

    # brca genes are concattenated
    if 'BRCA1' in concatted_symbols and 'BRCA2' in concatted_symbols:
        concatted_symbols.remove('BRCA1')
        concatted_symbols.remove('BRCA2')
        concatted_symbols.append('BRCA12')

    return concatted_symbols


def concatenate_files_with_identical_headers(file_directory, output):
    allFiles = glob.glob(file_directory + "/*.txt")
    with open(output, 'wb') as outfile:
        for i, fname in enumerate(allFiles):
            with open(fname, 'rb') as infile:
                if i != 0:
                    infile.readline()  # Throw away header on all but first file
                # Block copy rest of file from input to output without parsing
                shutil.copyfileobj(infile, outfile)


def get_lovd_symbols(symbols):
    lovd_symbols = []

    for symbol in concatenate_symbols(symbols):
        symbol = symbol.replace('BRCA12', 'BRCA')
        if symbol in ['BRCA', 'CDH1']:
            lovd_symbols.append(symbol)

    return lovd_symbols


def check_input_and_output_tsvs_for_same_number_variants(tsvIn, tsvOut,
                                                         numVariantsRemoved=0):
    tsvInput = csv.DictReader(open(tsvIn, 'r'), delimiter='\t')
    numVariantsIn = len(list(tsvInput))
    tsvOutput = csv.DictReader(open(tsvOut, 'r'), delimiter='\t')
    numVariantsOut = len(list(tsvOutput))
    print(
                "Number of variants in input: %s \nNumber of variants in output: %s \n Number of variants removed: %s\n" % (
        numVariantsIn, numVariantsOut, numVariantsRemoved))
    handle_process_success_or_failure(
        numVariantsIn - numVariantsRemoved == numVariantsOut, tsvOut)


def handle_process_success_or_failure(process_succeeded, file_path):
    file_name = file_path.split('/')[-1]
    if process_succeeded is True:
        print("Completed writing %s. \n" % (file_name))
    else:
        now = str(datetime.datetime.utcnow())
        file_directory = os.path.dirname(file_path)
        failed_file_name = "FAILED_" + now + "_" + file_name
        os.rename(file_path, file_directory + "/" + failed_file_name)
        print("**** Failure creating %s ****\n" % (file_name))


def extract_file(archive_path, tmp_dir, file_path):
    with tarfile.open(archive_path, "r:gz") as tar:
        tar.extract(file_path, tmp_dir)

    return tmp_dir + '/' + file_path

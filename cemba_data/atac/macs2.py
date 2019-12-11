import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def test_macs2():
    try:
        subprocess.run(['macs2', '--version'],
                       stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8',
                       check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print('Seems macs2 not correctly installed. Try to run: conda install "macs2>=2.2.4"')
        raise
    return


test_macs2()


def process_runner(cmd):
    p = subprocess.run(cmd,
                       stderr=subprocess.PIPE,
                       stdout=subprocess.PIPE,
                       encoding='utf8',
                       check=True)
    return p


def macs2(frag_bed_path_list, cpu, species, remove_temp=True, **macs2_kws):
    test_macs2()

    effective_genome_size = {
        'hs': '2.7e9',
        'mm': '1.87e9'
    }
    default_macs2_params = {
        '-g': effective_genome_size[species.lower()],
        "--nomodel": '',
        "--qval": "5e-2",
        "-B": '',
        "--SPMR": '',
        "--call-summits": '',
        "--keep-dup": "all",
        "--shift": "100",
        "--ext": "200",
        "-f": "AUTO"}
    default_macs2_params.update(macs2_kws)

    with ProcessPoolExecutor(cpu) as executor:
        futures = {}
        for path in frag_bed_path_list:
            path = str(path)
            output_prefix = path[:-7]  # remove .bed.gz

            macs2_params = default_macs2_params.copy()

            macs2_params['-t'] = path
            macs2_params['-n'] = output_prefix

            cmd = ['macs2', 'callpeak']
            for k, v in macs2_params.items():
                cmd.append(k)
                if v != '':
                    cmd.append(v)
            future = executor.submit(process_runner, cmd)
            futures[future] = output_prefix

        for future in as_completed(futures):
            try:
                future.result()
                output_prefix = futures[future]
                if remove_temp:
                    subprocess.run(['rm', '-f',
                                    f'{output_prefix}_control_lambda.bdg',
                                    f'{output_prefix}_peaks.xls',
                                    f'{output_prefix}_summits.bed',
                                    f'{output_prefix}_treat_pileup.bdg'])
            except subprocess.CalledProcessError as e:
                print(e.stderr)
                raise e
    return

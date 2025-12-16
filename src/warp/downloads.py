import paramiko
import os
import pandas as pd
import glob
from warp.server_config import cor_ips
from getpass import getpass


def build_server_path(raw_file, file_type='CCF'):
    split_f = raw_file.split('/')
    drs = split_f[1]
    instrument = split_f[0].upper() if '3.0.1' in drs else 'CORALIE14'
    drs = drs+'-DEV' if '3.0.1' in drs else drs
    path = f"/data/{instrument}DRS/{drs}/{split_f[2]}/{split_f[3]}/{split_f[4]}"
    if file_type.lower() == 'ccf':
        if '3.0.1' in drs:
            path = path.replace('.fits', '_CCF_A.fits')
        else:
            path = path.replace('.fits', '_ccf_*_A.fits')
    elif file_type.lower() == 's1d':
        if '3.0.1' in drs:
            path = path.replace('.fits', '_S1D_A.fits')
        else:
            path = path.replace('.fits', '_s1d_A.fits')
    elif file_type.lower() == 's2d':
        if '3.0.1' in drs:
            path = path.replace('.fits', '_S2D_BLAZE_A.fits')
        else:
            path = path.replace('.fits', '_e2ds_A.fits')
    else:
        raise ValueError(f"File type {file_type} not recognized.")
    return path


def download_files(file_list=None, file_type='CCF', output_directory='data/', verbose=True, password=None, user=None):
    """
    Download specified files from remote server via SSH tunneling.

    Args:
        drs (str): DRS version ('3.0.1' or '3.8').
        instrument (str): Instrument name.
    """
    if file_list is None:
        raise ValueError("file_list must be provided.")
    # We need to create a separate ssh connection for each coralie machine
    raw_list = {
        'coralie98': [],
        'coralie07': [],
        'coralie14': []
    }
    for file in file_list:
        raw_list[file.split('/')[0]].append(
            build_server_path(file, file_type=file_type))
    for instrument in raw_list.keys():
        if len(raw_list[instrument]) == 0:
            continue
        download_instrument_files(
            instrument, file_list=raw_list[instrument], save_dir=output_directory, verbose=verbose, user=user, password=password)
    return


def download_instrument_files(instrument, file_list=None, save_dir='data/', verbose=True, user=None, password=None):
    import fnmatch

    if instrument not in cor_ips.keys():
        raise ValueError(f"Instrument {instrument} not recognized.")
    if file_list is None or len(file_list) == 0:
        raise ValueError("file_list must be provided and non-empty.")
    host_ip = cor_ips[instrument]
    jump_host = "login01.astro.unige.ch"
    if user is None:
        user = input('Enter your username : ')
    if password is None:
        password = getpass('Enter your password : ')
    jump_client, final_client, sftp = open_connection(
        jump_host, host_ip, user, password, password, verbose=verbose)
    for remote_path in file_list:
        if '3.0.1' not in remote_path:
            # handle wildcards
            dir_path = os.path.dirname(remote_path)
            pattern = os.path.basename(remote_path)
            try:
                files = sftp.listdir(dir_path)
            except IOError:
                print(f"[WARN] Directory {dir_path} does not exist on server.")
                continue
            matched_files = fnmatch.filter(files, pattern)
            if len(matched_files) == 0:
                print(
                    f"[WARN] No files matching {pattern} found in {dir_path}.")
                continue
            remote_path = remote_path.replace(pattern, matched_files[0])
        basename = os.path.basename(remote_path)
        local_path = os.path.join(save_dir, instrument,
                                  remote_path.split('/')[-2],
                                  remote_path.split('/')[-1])
        if not os.path.exists(os.path.dirname(local_path)):
            os.makedirs(os.path.dirname(local_path), exist_ok=True)
        if (os.path.exists(local_path)):
            print(f"File {local_path} already exists. Skipping download.")
            continue
        if verbose:
            print(f"Downloading {remote_path} → {local_path}")
        sftp.get(remote_path, local_path)
        if verbose:
            print(f"✅ Downloaded {basename} to {local_path}")

    sftp.close()
    final_client.close()
    jump_client.close()


def open_connection(jump_host, final_host, user, jump_password, final_password, verbose=True):
    if verbose:
        print(f"[INFO] Connecting to jump host: {jump_host}")
    jump_client = paramiko.SSHClient()
    jump_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    jump_client.connect(jump_host, username=user, password=jump_password)

    jump_transport = jump_client.get_transport()
    dest_addr = (final_host, 22)
    local_addr = (jump_host, 22)
    channel = jump_transport.open_channel(
        "direct-tcpip", dest_addr, local_addr)
    if verbose:
        print(
            f"[INFO] Connecting to final host: {final_host} through jump host.")
    final_client = paramiko.SSHClient()
    final_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    final_client.connect(final_host, username=user,
                         password=final_password, sock=channel)

    # --- Open SFTP session to B ---
    sftp = final_client.open_sftp()

    return jump_client, final_client, sftp

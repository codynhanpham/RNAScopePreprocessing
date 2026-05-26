import os, sys, json
import hashlib

import typing


HASHFILENAME = "_datahashes.json"

BUF_CHUNK_SIZE = 1024 * 1024 * 10

data_default_schema = {
    "files": {},
    "vars": {}
}


def getHashFilePath(experiment_output_dir: str, experiment_name: str) -> str:
    return os.path.join(experiment_output_dir, experiment_name + HASHFILENAME)


def hashFile(filepath: str, chunk_size: int = BUF_CHUNK_SIZE) -> str:
    """
    Returns the sha256 hash of a file given its path.

    Parameters
    ----------
    filepath : str
        Path to the file.
    chunk_size : int, optional
        Size of the chunks to read the file, by default `BUF_CHUNK_SIZE`.

    Returns
    ----------
    `str`
        sha256 hash of the file
    """
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while True:
            data = f.read(chunk_size)
            if not data:
                break
            sha256.update(data)
    return sha256.hexdigest()


def hashVar(var_value: any) -> str:
    """
    Returns the sha256 hash of a variable.

    Parameters
    ----------
    var_value : any
        Value of the variable.

    Returns
    ----------
    `str`
        sha256 hash of the variable
    """
    sha256 = hashlib.sha256()
    sha256.update(str(var_value).encode())
    return sha256.hexdigest()



def updateHashFileLog(hashfilepath: str, hashkey: str, hash: str, type: typing.Literal["files", "vars"]) -> None:
    """
    Update the hash file with a new hash or add a new hash entry.

    Parameters
    ----------
    hashfilepath : str
        Path to the hash file.
    hashkey : str
        Key of the hash entry. Either the relative file path or the literal variable name.
    hash : str
        Hash to be added to the hash file.
    type : str
        Type of the hash key to look up. Either `"files"` or `"vars"`.
    """

    if os.path.exists(hashfilepath):
        with open(hashfilepath, 'r') as f:
            data = json.load(f)
    else:
        data = data_default_schema
        
    if type not in data:
        data[type] = {}
    data[type][hashkey] = hash

    with open(hashfilepath, 'w') as f:
        json.dump(data, f, indent=4)

    return None


def checkHashExists(hashfilepath: str, lookup_key: str, hash: str, type: typing.Literal["files", "vars"]) -> bool:
    """
    Check if a hash exists in a hash file, and that it is equal to the hash provided.
    
    Parameters
    ----------
    hashfilepath : str
        Path to the hash file.
    lookup_key : str
        Key to look up in the hash file.
    hash : str
        Hash to compare against the hash file.
    type : str
        Type of the hash key to look up. Either `"files"` or `"vars"`.

    Returns
    ----------
    bool
        True if the hash exists and equal in the hash file, False otherwise.
    """

    if not os.path.exists(hashfilepath):
        return False
    with open(hashfilepath, 'r') as f:
        data = json.load(f)

    if type not in data:
        return False
    
    if lookup_key not in data[type]:
        return False
    
    return data[type][lookup_key] == hash



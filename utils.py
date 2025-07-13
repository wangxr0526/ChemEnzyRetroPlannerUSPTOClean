import csv
import signal

from rdkit import Chem


class timeout:
    """
    Function for counting time. If a process takes too long to finish, it will be exited by this function.
    """

    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds

        self.error_message = error_message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        signal.alarm(0)


def canonicalize_smiles(smi, clear_map=True):
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        if clear_map:
            [atom.ClearProp('molAtomMapNumber') for atom in mol.GetAtoms() if atom.HasProp('molAtomMapNumber')]
        return Chem.MolToSmiles(mol)
    else:
        return ''


def get_writer(output_name, header):
    # output_name = os.path.join(cmd_args.save_dir, fname)
    fout = open(output_name, 'w')
    writer = csv.writer(fout)
    writer.writerow(header)
    return fout, writer

import os
import logging

from ngs import NGS
from ngs.Read import Read


class DirectSRAExtractor:
    def dump_to_fasta(self, sra_file_path):
        handles = []

        sra_base = os.path.basename(sra_file_path).replace('.sra','')

        forward_name = '{}_1.fasta'.format(sra_base)
        reverse_name = '{}_2.fasta'.format(sra_base)
        unpaired_name = '{}.fasta'.format(sra_base)

        handles.append(open(unpaired_name,'w'))
        handles.append(open(forward_name,'w'))
        handles.append(open(reverse_name,'w'))
        num_paired = 0
        num_unpaired = 0

        logging.info("Starting FASTA conversion ..")

        last_read_id = None
        with NGS.openReadCollection(os.path.abspath(sra_file_path)) as run:
            logging.info("Opened")
            with run.getReadRange(1, run.getReadCount(), Read.all) as it:
                while it.nextRead():
                    while it.nextFragment():
                        name = it.getReadId()
                        if it.isPaired():
                            # Not sure how else to know what is fwd and what is
                            # rev
                            num_paired += 1
                            if name == last_read_id:
                                h = handles[2]
                            else:
                                h = handles[1]
                            last_read_id = name
                        else:
                            num_unpaired += 1
                            h = handles[0]
                        h.write(">{}\n{}\n".format(
                            name, it.getReadBases()
                        ))
        to_return = []
        for h in handles:
            h.close()
        if num_paired == 0:
            os.remove(forward_name)
            os.remove(reverse_name)
        else:
            to_return.append(forward_name)
            to_return.append(reverse_name)
        if num_unpaired == 0:
            os.remove(unpaired_name)
        else:
            to_return.append(unpaired_name)

        return to_return

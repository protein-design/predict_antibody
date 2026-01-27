import numpy as np
import pandas as pd
import sys
import os

src_dir = os.path.dirname(os.getcwd())
bioomics_dir = '/home/yuan/bio/bio_omics/src'
for _dir in (src_dir, bioomics_dir):
    if _dir not in sys.path:
        sys.path.append(_dir)

from bioomics import QueryComplex

class LoadData:

    @staticmethod
    def antibody():
        query = "select * from view_antibody;"
        df = QueryComplex(True).list_data(query, True)
        print('pdb', len(df['pdb_id'].unique()))
        print('chains', len(df['chain_id'].unique()))
        return df
    
    @staticmethod
    def imgt_regions():
        query = f"""
            select V.pdb_id, V.chain_id, A.chain_no,
                V.region_name, V.seq_from, V.seq_to, V.seq
            from align_vfrag V
            left join view_antibody A on V.chain_id = A.chain_id
            where A.model_no = 0
        ;"""
        df = QueryComplex(True).list_data(query, True)
        df['seq_len'] = df['seq_to'] - df['seq_from'] + 1
        print('pdb', df['pdb_id'].nunique())
        print('chains', df['chain_id'].nunique())

        vg = df.groupby(['pdb_id', 'chain_no', 'region_name'])
        vg = dict(tuple(vg))
        return df, vg
    
    @staticmethod
    def distance(table_name):
        query = f"""
            select D.pdb_id, D.combo_id, C.chain_combo, D.relative_pkl
            from {table_name} D
            left join meta_combo2 C on D.combo_id = C.combo_id
            where D.pdb_id in (
                select distinct pdb_id from view_antibody
            )
            AND D.relative_pkl is not null
        ;"""
        rows = QueryComplex(True).list_data(query)
        df = pd.DataFrame(rows)
        print('pdb', df['pdb_id'].nunique())
        print('combo', df['combo_id'].nunique())
        print('regions', len(rows))
        return rows

    @staticmethod
    def cal_overlapping(ref_bound:tuple, test_bound:tuple, offset:int):
        c = 0
        ref_start, ref_end = ref_bound
        test_start, test_end = test_bound
        for i in range(ref_start-offset, ref_end+offset):
            if test_start <= i < test_end:
                c += 1
        r = c/(ref_end - ref_start)
        return 1 if r > 1 else r

    @staticmethod
    def scan_regions(region_name, motifs, vg, offset:int=0):
        res = []
        # motifs determined by experimental data
        mg = motifs.groupby(['pdb_id','chain_no'])
        for (pdb_id, chain_no), sub in mg:
            key = (pdb_id, chain_no, region_name)
            if key in vg:
                region = vg[key].iloc[0]
                aln_start, aln_end = int(region['seq_from'])-1, int(region['seq_to'])
                overlap = []
                for i, row in sub.iterrows():
                    exp_start, exp_end = row['start'], row['end']
                    _overlap = LoadData.cal_overlapping(
                        (aln_start, aln_end),
                        (exp_start, exp_end),
                        offset
                    )
                    overlap.append(_overlap)
                sub['overlap'] = overlap
                sub = sub[sub['overlap']>0].sort_values('overlap', ascending=False)
                aln = region.to_dict()
                if len(sub):
                    exp = sub.iloc[0].to_dict()
                    aln.update({
                        'aa': exp.get('seq'),
                        'pair_aa': exp.get('pair_aa'),
                        'exp_start': exp.get('start'),
                        'exp_end': exp.get('end'),
                        'pair_chain_no': exp.get('pair_chain_no'),
                        'combo_id': exp.get('combo_id'),
                        'overlap': exp.get('overlap'),
                    })
                res.append(aln)
        print(region_name, len(res))
        return pd.DataFrame(res).fillna(0)
    
    @staticmethod
    def plddt():
        adf = LoadData.antibody()
        query = "select * from chain_plddt p;"
        pdf = QueryComplex(True).list_data(query, True)
        df = pd.merge(pdf, adf, how='left', on='chain_id')
        return df

    @staticmethod
    def plddt_rmsd():
        adf = LoadData.antibody()
        query = """
            select r.chain_id, r.rmsd, p.avg_plddt
            from chain_rmsd       r
            left join chain_plddt p
            on r.chain_id=p.chain_id and r.model_pdb = p.relative_pdb
            where r.chain_status = 'Raw'
                and r.rmsd is not null
                and p.ranking=1
        ;"""
        pdf = QueryComplex(True).list_data(query, True)
        df = pd.merge(pdf, adf, how='left', on='chain_id')
        return df

    @staticmethod
    def plddt_tm():
        adf = LoadData.antibody()
        query = """
            select r.chain_id, r.tm1, r.rmsd, p.avg_plddt
            from chain_tmalign r
            left join chain_plddt p
            on r.chain_id=p.chain_id and r.model_pdb = p.relative_pdb
            where r.chain_status = 'Raw'
                and r.tm1 is not null
                and p.ranking=1
        ;"""
        pdf = QueryComplex(True).list_data(query, True)
        df = pd.merge(pdf, adf, how='left', on='chain_id')
        return df
    
    @staticmethod
    def equal_seq():
        query = """
            select * from pdb_chain_seq
            where chain_id in (
                select chain_id from view_antibody
            )
        """
        df = QueryComplex(True).list_data(query, True)
        g = df.groupby('first_chain_id')
        return g
        
    @staticmethod
    def equal_seq_rmsd():
        adf = LoadData.antibody()
        query = """
            select * from equal_seq_rmsd
            where chain_id in (
                select chain_id from view_antibody
            )
        """
        df = QueryComplex(True).list_data(query, True)
        df = pd.merge(df, adf, how='left', on='chain_id')
        df = df.dropna()
        return df
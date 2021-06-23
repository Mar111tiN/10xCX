import os
import pandas as pd
import numpy as np
from script_utils_10x import show_output, cmd2df


def bam2df(bam_file, region="", bed="", config={}):
    """
    wrapper building CLI chain around the core tool bam2csv.mawk

    """

    # PARAMS
    # unwrap mawk tools
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    c = config["10x"]
    p = config["paths"]
    bam_path = os.path.join(os.getenv("HOME"), p["bampath"])
    # get the full bam file path
    if not bam_file.startswith("/"):
        bam_file = os.path.join(bam_path, bam_file)

    # create the basic command and unpack required params

    # ##### FILTERBED
    sam_cmd = f"samtools view {bam_file}"
    if region:
        sam_cmd += f" -r {region}"
    if bed:
        sam_cmd += f" -l {bed}"

    b2c_cmd = f"{mawk('bam2csv')} -c CHR,START,FLAG,MAPQ,CIGAR,LENGTH,SEQ -t CB,UB"

    cmd = f"{sam_cmd} | {b2c_cmd}"

    try:
        bam_df = cmd2df(cmd, show=True, multi=False)
    except Exception as e:
        show_output(f"There was an error using shell command <<{e}>>", color="warning")
        return cmd
    bam_df = bam_df.rename(
        dict(
            CHR="Chr",
            START="Start",
            FLAG="Flag",
            CIGAR="Cigar",
            LENGTH="SeqLen",
            SEQ="Seq",
        ),
        axis=1,
    )
    bam_df.loc[:, "SeqLen"] = bam_df["Seq"].str.len()
    # strip terminal "-1" in CB
    bam_df.loc[:, "CB"] = bam_df["CB"].str.replace("-1$", "", regex=True)

    # filter
    qmin = c["MAPQmin"]
    bam_df = bam_df.sort_values(["Start", "CB", "UB"]).query("MAPQ > @qmin")
    return bam_df


def get_intron_pos(df, intron=[0, 1]):
    """
    detect the read position relative to the intron:

    #### IntronPosAbs:
    0: Left
    1: spanningLeft
    2: within intron
    3: over intron
    4: spanningRight
    5: Right
    #### IntronPos:
    0: outside
    1: spanning exon-intron border
    2: within intron
    3: over int
    """

    cols = ["M1", "N", "M2"]
    # extract End Position from cigar string
    df.loc[:, cols] = df["Cigar"].str.extract(
        "^(?:[0-9]+S)?(?:(?P<M1>[0-9]+)M)(?:(?P<N>[0-9]+)N)?(?:(?P<M2>[0-9]+)M)?$"
    )
    for col in cols:
        df.loc[:, col] = df[col].fillna(0).astype(int)

    # write less do more
    i0 = intron[0]
    i1 = intron[1]
    # get absolute intron pos as accumulated intron bool conditions of Start and End coords
    df.loc[:, "GenLen"] = df["M1"] + df["N"] + df["M2"]
    df.loc[:, "End"] = df["Start"] + df["GenLen"]
    df.loc[:, "IntronPosAbs"] = (
        (df["Start"] > i0).astype(int)
        + (df["Start"] > i1).astype(int)
        + (df["End"] > i0).astype(int)
        + (df["End"] > i1).astype(int) * 2
    )
    df.loc[df["IntronPosAbs"] > 3, "IntronPos"] = 5 - df["IntronPosAbs"]
    df.loc[df["IntronPosAbs"] < 4, "IntronPos"] = df["IntronPosAbs"]
    df.loc[:, "IntronPos"] = df["IntronPos"].astype(int)
    base_cols = [
        "CB",
        "UB",
        "IntronPos",
        "IntronPosAbs",
        "Chr",
        "Start",
        "End",
        "Cigar",
        "GenLen",
        "Seq",
        "SeqLen",
        "MAPQ",
        "Flag",
    ]
    return df.loc[:, base_cols]


def get_Alt(df):
    """
    aggregate utility for row-wise accumulation of IntronPos info
    """
    s = pd.Series(
        dict(
            Alt=(df["IntronPos"] == 3).sum(),
            AB=((df["IntronPos"] > 0) & (df["IntronPos"] < 3)).sum(),
            Irr=(df["IntronPos"] == 0).sum(),
        )
    )
    return s


def detect_CXCR3Alt(bam_file, region="", bed_file="", intron=[], config={}):
    """
    master func
    """
    show_output(f"Loading bam file {bam_file}.")
    df = bam2df(bam_file, region=region, bed=bed_file, config=config)
    show_output(
        f"Detected a total of {len(df['UB'].unique())} molecules in {len(df['CB'].unique())} cells"
    )

    show_output(f"Looking for intron events in {len(df.index)} reads.")

    intron_df = get_intron_pos(df, intron=intron)
    show_output(
        f"Detected {len(intron_df.query('IntronPos > 0').index)} reads with intron involvement"
    )

    alt_df = intron_df.groupby(["CB", "UB"]).apply(get_Alt).reset_index()
    t_df = alt_df.query("Alt > 0 or AB > 0")
    show_output(
        f"Detected a total of {len(t_df['UB'].unique())} molecules with intron involvement in {len(t_df['CB'].unique())} cells",
        color="success",
    )
    return t_df

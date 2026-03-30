## This is all be tested

import altair as alt
import polars as pl
import yaml

# if analysis == "heterozygosity":
#     input_file = config["results"] + 'heterozygosity/gvcf_counts.het'
#     indiv_col = "SAMPLE"
#     value_cols = ["HET"]

# if analysis == "inbreeding":
#     input_file = config["results"] + f"inbreeding/GCTA/pass/{config['dataset']}.{config['subset']}.ibc"
#     indiv_col = "IID"
#     value_cols = ["NOMISS", "Fhat1", "Fhat2", "Fhat3"]

    # # Output path
    # output_dir = "/master/abagwell/variant-analysis/results/rhesus/heterozygosity/plots"

# Read config file. User should update this file path as needed
configfile = "/master/abagwell/workspace/github_project/variant-analysis/config/rhesus.yaml"
with open(configfile, 'r') as file:
    config = yaml.safe_load(file)

    # Load colors
    colors = pl.read_csv(config["colors"], separator="\t", schema_overrides={"Cohort": pl.String})

    # Load cohorts
    colonies_file = config["cohorts"]
    colonies = pl.read_csv(colonies_file, separator="\t", comment_prefix="#", schema_overrides={"Id": pl.String}).select("Id", "Cohort").unique()

def load_tsv(input_file, indiv_col, value_cols, seq_type="WGS"):
    """Load file and do initial processing.
    
    seq_type: Seq type to filter on. WES and WGS are not compatible. WES has higher heterozygosity
    """

    return pl.read_csv(input_file,
        separator='\t',
        schema_overrides={indiv_col: pl.String}
    ).with_columns(
        Seq = pl.col(indiv_col).str.slice(0, 3),
        #   ).str.strip_prefix("WES"
        #   ).str.strip_prefix("WGS"
        #   ).str.strip_prefix("AMP"
        #   ).str.strip_prefix("GBS"
        Indiv = pl.col(indiv_col).str.split('_').list.get(0).str.slice(3),
    ).drop(indiv_col
    # Filter to one type of sequencing
    ).filter(
        pl.col("Seq") == seq_type
    # Deduplicate individuals
    ).group_by('Indiv').agg(pl.first('*')
    # Join cohort info
    ).join(colonies, left_on='Indiv', right_on='Id'
    ).with_columns(
        pl.col('Cohort').cast(pl.Enum( 
            list(colors["Cohort"])
        )),
    ).group_by("Cohort").agg(*value_cols, 'Indiv'
    ).with_row_index("pop_idx", offset=1).with_columns(
        # Find which are year ranges
        is_year = pl.col("Cohort").cast(pl.String).str.contains("-").not_().cast(pl.Int8)
    ).with_columns(
        # Set the index of year ranges to 0
        pl.col("pop_idx").mul("is_year")
    ).drop("is_year").sort("Cohort"
    ).group_by(
        # Create color index
        "pop_idx", maintain_order=True
    ).agg('*').with_row_index("color_idx"
    ).drop("pop_idx"
    ).explode(pl.exclude("color_idx")
    ).explode(pl.exclude("Cohort", "color_idx"))


def df_with_trios(df, value_col):
    # Read pedigree to find trios
    pedigree = pl.read_csv(config["demographics"], separator="\t", columns=['Id', 'Sire', 'Dam'],
        schema_overrides={'Id': pl.String, 'Sire': pl.String, 'Dam': pl.String})

    # Create dataframe for plotting trios
    return df.join(
        # Add pedigree information
        pedigree, left_on='Indiv', right_on='Id', how='left'
    ).filter(
        # Keep only P51 offspring
        pl.col("Cohort") == "P51 offspring"  # TODO: Generalize
    ).sort(value_col, descending=True).with_row_index("trio_idx"
    ).unpivot(["Indiv", "Sire", "Dam"], index="trio_idx", variable_name="Relation", value_name="Id").join(
        # Rejoin
        df, left_on='Id', right_on='Indiv', how='left'
    )

def add_parental_cohorts(df, indiv_col, value_col, original_df):
    # Find U42 sires
    sires_df = df.filter(pl.col("Relation") == "Sire").select("Id").sort("Id").with_columns(
        Cohort = pl.lit("U42 sires")
    )

    # Find CPRC dams
    dams_df = df.filter(pl.col("Relation") == "Dam").select("Id").sort("Id").with_columns(
        Cohort = pl.lit("CPRC dams")
    )

    return pl.concat([colonies, sires_df, dams_df]
    ).join(original_df.select("color_idx", indiv_col, value_col), left_on="Id", right_on=indiv_col, how="left"
    ).filter(pl.col(value_col).is_not_null()).unique().rename({"Id": "Indiv"})

def plot_boxplot(df, value_col):
    return alt.Chart(df).mark_boxplot().encode(
        alt.X('Cohort', title='Cohort', sort=colors["Cohort"]),
        alt.Y(value_col, title=value_col).scale(zero=False), 
        alt.Color('Cohort:N', legend=None,).scale(
            domain = list(colors["Cohort"]),
            range = list(colors["Color"])
        ),
        tooltip=[
            alt.Tooltip('Indiv')
        ]
    ).properties(
        title=value_col
    )

def plot_violinplot(df, df_with_parental_cohorts, value_col):
    """Plot narrow violinplot."""

    # Varibles to adjust
    error_unit = 'stderr' # Can switch extent to `stdev`, `stderr`, or `ci`

    # Find ceiling and floor for plot
    max_y = df.select(pl.col(value_col).max().mul(1000).ceil().truediv(1000)).item()
    min_y = df.select(pl.col(value_col).min().mul(1000).floor().truediv(1000)).item()

    violin = alt.Chart().transform_density(
        value_col,
        as_=[value_col, 'density'],
        extent=[min_y, max_y],
        groupby=['Cohort']
    ).mark_area(orient='horizontal').encode(
        alt.X("density:Q").stack('center')
            .impute(None)
            .title(None)
            .axis(labels=False, values=[0], grid=False, ticks=True)
            .scale(nice=False,zero=False),
        alt.Y(f"{value_col}:Q", title=value_col),#.scale(domain=[min_y, max_y]),
        # alt.Column("Cohort:N", title="Cohort",
        #       # TODO: Generalize this
        #     ).spacing(0).header(titleOrient='bottom', labelOrient='bottom', labelPadding=0),
        #color=alt.Color("color_idx:N", legend=None).scale(scheme="category10"),
        color=alt.Color("Cohort:N", legend=None).scale(
            domain = list(colors["Cohort"]),
            range = list(colors["Color"])
        )
    ).properties(
        width=25
    )

    error = alt.Chart().mark_errorbar(extent=error_unit).encode(
        #alt.X('Cohort', title=None),
        alt.Y(f'{value_col}:Q', title=value_col)
        )

    mean = alt.Chart().mark_circle(color='black').encode(
        #alt.X('Cohort', title=None),
        alt.Y(f'mean({value_col}):Q', title=value_col)
        )

    return alt.layer(violin, error, mean, data=df_with_parental_cohorts
        # .filter(
        # #     # pl.col("Cohort").is_in(["2018-2020", "Offspring of merger", "NEPRC source"]))
        #     pl.col("Cohort").is_in(["Conventional source", "Brooks source", "NEPRC source"]))
    ).facet(
        #column='Cohort'
        alt.Column('Cohort',
            header=alt.Header(
                labelOrient='bottom', labelPadding=0, labelAnchor='middle', labelAngle=-90, labelBaseline="middle", labelAlign="right",
                title='Cohort', titleAlign="center", titleOrient='bottom', ),
            sort=colors["Cohort"]
        )
    ).resolve_scale(x=alt.ResolveMode("independent")
    ).configure_facet(
        spacing=0,
    ).configure_title(anchor='middle').properties(
        title=[f"{value_col} in Cohorts"]
    )

def plot_lineplot(df, final_value_col, inverse=False):
    # Remove "U42 sires" since they are already included under "U42"
    simplified_colors = colors.filter(~pl.col("Cohort").is_in(["U42 sires", "CPRC dams"]))

    max_y = df.select(pl.col(final_value_col).max().mul(1000).ceil().truediv(1000)).item()
    min_y = df.select(pl.col(final_value_col).min().mul(1000).floor().truediv(1000)).item()

    # Plot trios to more easily compare offspring to parents
    return alt.Chart(df).mark_line().encode(
        alt.X("trio_idx", title='Trio',
              axis=alt.Axis(labels=False, ticks=False,),
              scale=alt.Scale(domainMax=df.select("trio_idx").max().item()),
              sort="descending" if inverse else "ascending"
              ),
        alt.Y(final_value_col, title=final_value_col, scale=alt.Scale(domain=[min_y, max_y])),
        color=alt.Color("Cohort:N").scale(
            domain = list(simplified_colors["Cohort"]),
            range = list(simplified_colors["Color"])
            ),
        tooltip=[
            alt.Tooltip('Id'),
            alt.Tooltip('Relation')
        ]
    ).properties(
        title=f"Comparison of {final_value_col} in P51 Offspring and Parents"
    )#.save(f"{output_dir}/P51_offspring_trios.heterozygosity.lineplot.html")
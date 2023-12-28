#' @title GENI plots test data for the Manhattan plot
#'
#' @description An R object containing formatted results from a genome-wide
#'   association study (GWAS) of interleukin-6 levels with *p*-values < `1e-3`
#'   ([GCST90274815](https://www.ebi.ac.uk/gwas/studies/GCST90274815)).
#'
#' @name geni_test_manhattan
#'
#' @format A `data.frame` containing the variables necessary for
#'   constructing a static Manhattan plot. Specifically, this includes
#'   the following information:
#'
#' \itemize{
#'   \item{\code{chr}} {
#'     the chromosome for each genomic marker
#'   }
#'   \item{\code{pos}} {
#'     the genomic position (GRCh37) for each genomic marker
#'   }
#'   \item{\code{pvalue}} {
#'     the association p-value for each genomic marker
#'   }
#'   \item{\code{highlight}} {
#'     an indicator variable whether the genomic marker should be highlighted:
#'     `0` = point not highlighted,
#'     `1` = first highlight colour
#'   }
#'   \item{\code{highlight_shape}} {
#'     an indicator variable of the illustrated shape for the genomic markers:
#'     `0` = standard circle,
#'     `5` = standard diamond with border
#'   }
#'   \item{\code{label}} {
#'     the name of the nearest gene for each sentinel genome-wide
#'     associated genetic variant, the labels for non-sentinel
#'     genetic variants are omitted
#'   }
#' }
#'
#' @examples
#' head(geni_test_manhattan)
#'
#' @author James Staley <jrstaley95@gmail.com>
#' @author Wes Spiller
#'
#' @md

"geni_test_manhattan"


#' @title GENI plots test data for the PheWAS plot
#'
#' @description An R object containing formatted data from a phenome-wide
#'   association study (PheWAS) of rs2228145
#'   ([1:154426970-A-C](https://r9.finngen.fi/variant/1:154454494-A-C))
#'   from FinnGen r9.
#'
#' @name geni_test_phewas
#'
#' @format A `data.frame` containing the variables necessary for
#'   constructing an interactive PheWAS plot. Specifically, this includes
#'   the following information:
#'
#' \itemize{
#'   \item{\code{pvalue}} {
#'     the association p-value for each phenotype
#'   }
#'   \item{\code{sign}} {
#'     the direction of the association with the phenotype, where
#'     `0` = missing,
#'     `1` positive association,
#'     `-1` negative association
#'   }
#'   \item{\code{group}} {
#'     the phenotype group for each phenotype
#'   }
#'   \item{\code{label}} {
#'     the phenotype name
#'   }
#'   \item{\code{text}} {
#'     the hover text containing further information on the association
#'     including: phenotype, phenotypic category, genetic variant ID,
#'     direction of association, association p-value, number of cases
#'     and number of controls
#'   }
#' }
#'
#' @examples
#' head(geni_test_phewas)
#'
#' @author James Staley <jrstaley95@gmail.com>
#' @author Wes Spiller
#'
#' @md

"geni_test_phewas"

#' @title GENI plots test data for the regional plot
#'
#' @description An R object containing formatted results from a genome-wide
#'   association study (GWAS) of interleukin-6 levels
#'   ([GCST90274815](https://www.ebi.ac.uk/gwas/studies/GCST90274815))
#'   and linkage disequilibrium statistics from
#'   [1000 Genomes Phase 3](https://www.internationalgenome.org/)
#'   for the region 1:154301970-154551970.
#'
#' @name geni_test_region
#'
#' @format A `list` containing 2 objects with the information necessary
#'   for constructing an interactive regional plot.
#'   Specifically, this includes the following information:
#'
#' \itemize{
#'   \item{\code{assoc}} {
#'     a `data.frame` with genetic association results for interleukin-6 levels
#'     from [GCST90274815](https://www.ebi.ac.uk/gwas/studies/GCST90274815) for
#'     the region 1:154301970-154551970
#'   }
#'   \item{\code{corr}} {
#'     a `matrix` of Pearson correlation statistics (not squared) from the
#'     European samples of
#'     [1000 Genomes Phase 3](https://www.internationalgenome.org/)
#'     for the region 1:154301970-154551970, this `matrix` has the same markers
#'     in the same order as `assoc`
#'   }
#' }
#'
#' @details
#' \code{assoc} is a `data.frame` with the following columns:
#' \itemize{
#'   \item{\code{marker}} {
#'     the genomic marker identifier (i.e. rsID or chromosome-position)
#'   }
#'   \item{\code{chr}} {
#'     the chromosome for each genomic marker
#'   }
#'   \item{\code{pos}} {
#'     the genomic position (GRCh37) for each genomic marker
#'   }
#'   \item{\code{pvalue}} {
#'     the association p-value for each genomic marker
#'   }
#' }
#'
#' @author James Staley <jrstaley95@gmail.com>
#' @author Wes Spiller
#'
#' @examples
#' head(geni_test_region$assoc)
#' head(geni_test_region$corr)
#'
#' @md

"geni_test_region"

#' @title GENI plots test data for the stacked regional plot
#'
#' @description An R object containing formatted results from genome-wide
#'   association studies (GWAS) of interleukin-6 levels
#'   ([GCST90274815](https://www.ebi.ac.uk/gwas/studies/GCST90274815)) and
#'   interleukin-6 receptor levels
#'   ([GCST90088597](https://www.ebi.ac.uk/gwas/studies/GCST90088597))
#'   and linkage disequilibrium statistics from
#'   [1000 Genomes Phase 3](https://www.internationalgenome.org/)
#'   for the region 1:154301970-154551970.
#'
#' @name geni_test_stack_region
#'
#' @format A `list` containing 2 objects with the information necessary
#'   for constructing an interactive stacked regional plot.
#'   Specifically, this includes the following information:
#'
#' \itemize{
#'   \item{\code{assoc}} {
#'     a `data.frame` with genetic association results for interleukin-6 levels
#'     from [GCST90274815](https://www.ebi.ac.uk/gwas/studies/GCST90274815) and
#'     for interleukin-6 receptor levels
#'     [GCST90088597](https://www.ebi.ac.uk/gwas/studies/GCST90088597)
#'     for the region 1:154301970-154551970
#'   }
#'   \item{\code{corr}} {
#'     a `matrix` of Pearson correlation statistics (not squared) from the
#'     European samples of
#'     [1000 Genomes Phase 3](https://www.internationalgenome.org/)
#'     for the region 1:154301970-154551970, this `matrix` has the same markers
#'     in the same order as `assoc`
#'   }
#' }
#'
#' @details
#' \code{assoc} is a `data.frame` with the following columns:
#' \itemize{
#'   \item{\code{marker}} {
#'     the genomic marker identifier (i.e. rsID or chromosome-position)
#'   }
#'   \item{\code{chr}} {
#'     the chromosome for each genomic marker
#'   }
#'   \item{\code{pos}} {
#'     the genomic position (GRCh37) for each genomic marker
#'   }
#'   \item{\code{pvalue_1}} {
#'     the association p-value for interleukin-6 levels for each genomic
#'     marker
#'   }
#'   \item{\code{pvalue_2}} {
#'     the association p-value for interleukin-6 receptor levels for each
#'     genomic marker
#'   }
#' }
#'
#' @author James Staley <jrstaley95@gmail.com>
#' @author Wes Spiller
#'
#' @examples
#' head(geni_test_stack_region$assoc)
#' head(geni_test_stack_region$corr)
#'
#' @md

"geni_test_stack_region"

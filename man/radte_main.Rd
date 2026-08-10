\name{radte_main}
\alias{radte_main}
\title{Run the RADTE command-line workflow}
\usage{
radte_main(command_args = commandArgs(trailingOnly = TRUE))
}
\arguments{
  \item{command_args}{A character vector of RADTE command-line arguments in
  \code{--key=value} form.}
}
\value{
Invisibly returns a list describing the output directory, prefix, manifest,
dating backend, and calibrated node class after a successful run. Help and
version requests invisibly return zero.
}
\description{
Validates command-line arguments, loads reconciled gene-tree and dated
species-tree inputs, runs the selected dating backend, and writes RADTE
artifacts plus an audit manifest.
}

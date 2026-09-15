# Frozen scalar-loop reference

These six files preserve the direct maker used before the RDF migration:
maker, reader, selection, weights, configuration and I/O. Configuration and I/O
were frozen at `b32c0eb` during the direct-maker simplification so new production
record types cannot change this oracle. Include paths and provenance file paths
were adjusted to this directory; numerical logic is unchanged.

Run the reference in a separate ROOT process from the current maker: both use
the `PJ` namespace. The reference supplies no-smearing comparisons; historical
RDF outputs and candidate audits supply the accepted keyed-smearing realization.
Its old job-order response RNG is deliberately not the production RNG.

Normal workflows use `histmakers/PhotonJetHistMaker.C`. Historical RDF and hybrid
benchmark implementations remain identifiable by their preserved Git revisions.

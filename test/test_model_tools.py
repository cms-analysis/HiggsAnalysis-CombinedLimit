#!/usr/bin/env python3
"""Check that ModelBuilderBase.factory_ parses RooWorkspace statements correctly."""

from types import SimpleNamespace

import ROOT

from HiggsAnalysis.CombinedLimit.ModelTools import ModelBuilderBase


for n_terms in (1, 400):
    # Setup
    builder = ModelBuilderBase.__new__(ModelBuilderBase)
    builder.options = SimpleNamespace(verbose=0)
    builder.out = ROOT.RooWorkspace("w", "w")
    builder.out.dont_delete = []

    # This sums n_terms copies of r and adds s.
    # With 400 terms the statement exceeds 1000 characters; the final s term
    # checks that parsing reaches the end of the formula.
    variable = builder.factory_("r[2,0,10]")
    tail_variable = builder.factory_("s[5,0,10]")
    formula = "+".join(["@0"] * n_terms) + "+@1"
    statement = f'expr::expression("{formula}",r,s)'
    if n_terms == 400:
        assert len(statement) > 1000
    expression = builder.factory_(statement)

    assert variable == builder.out.var("r")
    assert tail_variable == builder.out.var("s")
    assert expression
    assert expression == builder.out.function("expression")
    assert expression.getVal() == n_terms * 2.0 + 5.0, expression.getVal()
    # Changing both inputs also checks that the end of the statement was parsed.
    variable.setVal(3.0)
    tail_variable.setVal(7.0)
    assert expression.getVal() == n_terms * 3.0 + 7.0, expression.getVal()
    assert len(builder.out.dont_delete) == 3
    assert builder.out.dont_delete[0] is variable
    assert builder.out.dont_delete[1] is tail_variable
    assert builder.out.dont_delete[2] is expression

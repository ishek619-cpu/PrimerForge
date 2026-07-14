"""
Workflow tests.
"""

from envoprimer.species.workflow import WorkflowResult


def test_workflow_result():

    result = WorkflowResult(

        dataset=None,

        alignment=None,

        diagnostic_sites=[],

        diagnostic_windows=[],

    )

    assert result.dataset is None

    assert result.alignment is None

    assert result.diagnostic_sites == []

    assert result.diagnostic_windows == []

"""Search engine over diameter distributions for ``rcpgenerator``.

User-facing surface (re-exported here for one-import convenience):

* MODEL: :func:`make_universal_model`, :func:`balance_amplitudes`,
  :func:`set_perturb_weights`, :func:`initialize_theta`, :data:`COMPONENT_REGISTRY`
* Packing: :func:`create_packing`, :func:`run_one_custom_packing`
* Search: :func:`run_or_resume_search`, :func:`run_validation`
* State + accessors: :class:`SearchState`
* Diagnostics: :func:`run_all_diagnostics`, :func:`plot_search_history`,
  :func:`plot_distribution_from_theta`, :func:`plot_distribution_evolution`,
  :func:`plot_population_at_generation`, :func:`regenerate_packing`,
  :func:`compare_runs`
* Quality: :func:`overlap_report`, :func:`phi_corrected`
* Voronoi: :func:`voronoi_phi_local`, :func:`plot_local_phi_histogram`,
  :func:`plot_local_phi_vs_diameter`
* Studies: :func:`list_studies`, :func:`build_study`, :func:`default_config`
"""
from ._helpers import (
    auto_neighbor_max,
    get_final_phi,
    hypersphere_volume_from_diameter,
    safe_compute_phi,
    summarize_packing,
)
from .diagnostics import (
    compare_runs,
    plot_distribution_evolution,
    plot_distribution_from_theta,
    plot_population_at_generation,
    plot_search_history,
    regenerate_packing,
    run_all_diagnostics,
)
from .jobs import (
    aggregate_candidate_results,
    create_packing,
    make_packing_job,
    resolve_n_workers,
    run_jobs,
    run_one_custom_packing,
    run_packing_job,
)
from .model import (
    COMPONENT_REGISTRY,
    balance_amplitudes,
    center_for_geometric_mean,
    clip_theta,
    evaluate_model_pdf,
    initialize_theta,
    lognormal_diagnostics,
    make_universal_model,
    perturb_theta,
    sample_diameters,
    set_perturb_weights,
    theta_to_dataframe,
    tukey_window_on_interval,
    universal_pdf_R,
)
from .overlap import overlap_report, phi_corrected
from .search_v2 import run_or_resume_search_v2
from .snapshot import restore_from_snapshot, load_snapshot
from .census import load_census, moments_to_delta_S, census_path
from .sweep import (
    build_sweep_model,
    autosize_N,
    make_grid,
    plan_sweep,
    run_one_case,
    run_sweep_parallel,
    run_sweep_sequential,
    load_sweep_census,
)
# optional targeted save/bundle for the webapp + analysis (lazy deps inside)
from .bundling import bundle_from_npz, save_outputs, save_npz
from .persistence import (
    SearchState,
    load_config_sidecar,
    load_model_sidecar,
    save_config_sidecar,
    save_model_sidecar,
)
from .search import (
    propose_population,
    random_theta,
    run_one_generation,
    run_or_resume_search,
    run_validation,
)
from .studies import (
    STUDIES_2D,
    STUDIES_3D,
    build_study,
    default_config,
    list_studies,
)

# Voronoi is optional — only importable when pyvoro-mmalahe is installed.
try:
    from .voronoi import (
        plot_local_phi_histogram,
        plot_local_phi_vs_diameter,
        voronoi_phi_local,
    )
except ImportError:
    voronoi_phi_local = None
    plot_local_phi_histogram = None
    plot_local_phi_vs_diameter = None

__all__ = [
    # model
    "COMPONENT_REGISTRY", "make_universal_model", "balance_amplitudes",
    "set_perturb_weights", "theta_to_dataframe", "initialize_theta", "evaluate_model_pdf",
    "sample_diameters", "clip_theta", "perturb_theta", "universal_pdf_R",
    "tukey_window_on_interval", "lognormal_diagnostics", "center_for_geometric_mean",
    # bundling (optional save/bundle for webapp + analysis)
    "bundle_from_npz", "save_outputs", "save_npz",
    # jobs
    "create_packing", "run_one_custom_packing", "run_packing_job",
    "make_packing_job", "run_jobs", "resolve_n_workers",
    "aggregate_candidate_results",
    # search
    "run_or_resume_search", "run_one_generation", "run_validation",
    "propose_population", "random_theta",
    # persistence
    "SearchState", "save_config_sidecar", "save_model_sidecar",
    "load_config_sidecar", "load_model_sidecar",
    # diagnostics
    "run_all_diagnostics", "plot_search_history",
    "plot_distribution_from_theta", "plot_distribution_evolution",
    "plot_population_at_generation", "regenerate_packing", "compare_runs",
    # quality
    "overlap_report", "phi_corrected",
    "run_or_resume_search_v2",
    "restore_from_snapshot", "load_snapshot",
    "load_census", "moments_to_delta_S", "census_path",
    # continuous-baseline sweeps
    "build_sweep_model", "autosize_N", "make_grid", "plan_sweep",
    "run_one_case", "run_sweep_parallel", "run_sweep_sequential",
    "load_sweep_census",
    # voronoi
    "voronoi_phi_local", "plot_local_phi_histogram",
    "plot_local_phi_vs_diameter",
    # studies
    "STUDIES_2D", "STUDIES_3D", "list_studies", "build_study",
    "default_config",
    # helpers
    "auto_neighbor_max", "get_final_phi", "safe_compute_phi",
    "hypersphere_volume_from_diameter", "summarize_packing",
]

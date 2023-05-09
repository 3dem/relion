import typer

CLI_NAME = 'relion_tomo_derive_particle_positions'
cli = typer.Typer(
    name=CLI_NAME,
    add_completion=False,
    no_args_is_help=True,
    help='Derive particle poses from annotations for refinement.'
)

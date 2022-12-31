$simulators = "first_order_upwind", "lax_wendroff",  "beam_warming", "fromm", "tvd_minmod", "tvd_superbee", "tvd_van_leer", "tvd_van_albada"
foreach ($simulator in $simulators)
{
    & ./build/$simulator
}
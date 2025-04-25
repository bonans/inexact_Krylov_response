include("silicon.jl")
debug = false
setup_model(debug ? "debug" : "full")
run_gmres(; debug=debug, restart=20, tol=1e-9, adaptive=true, CG_tol_scale_choice="agr", precon=false)
run_gmres(; debug=debug, restart=20, tol=1e-9, adaptive="D10", CG_tol_scale_choice=nothing, precon=false)
run_gmres(; debug=debug, restart=20, tol=1e-9, adaptive=true, CG_tol_scale_choice="hdmd", precon=false)
run_gmres(; debug=debug, restart=20, tol=1e-9, adaptive="D100", CG_tol_scale_choice=nothing, precon=false)
run_gmres(; debug=debug, restart=20, tol=1e-9, adaptive=true, CG_tol_scale_choice="grt", precon=false)
run_gmres(; debug=debug, restart=20, tol=1e-9, adaptive="D10_n", CG_tol_scale_choice=nothing, precon=false)
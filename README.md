# Nates-OSPM
Package Layout
Gal_Dynamics
	Data
		Data_Prep 
			DATA_NEW_GAL.py
			Data_Asmbl.py
			Data_Config.py (NEW NEEDS WORK)
			Data_Geometry.py
			Data_Load_Center.py
			Data_Paths.py
			Data_Preprocess.py
			Data_Sources.py
			Data_Web_Srcs.py
			Data__README.md
			build_centers.sh
			Galaxies.txt 
		Profiles 
			All the Galaxy Data
				Dump
				Profiles
					Galaxy_Stellar_Config.py
					Galaxy_OSPM_Config.py	
				checkpoint
				data
				logs
			Draco
				Profiles
					Draco_Stellar_Profile.csv
					Draco_OSPM_Config.py
				Dump
				Testing (archive later)
				checkpoint
				data
					Draco.csv (I think this one is the same as Draco_Stellar_Profile.csv)
					README.txt
				logs 
		center(Merging w/ Profiles)
	OSPM
		AI
			OSPM_Daemon.py
		Controllers
			OSPM_API.py
			OSPM_Control.py
			OSPM_MASTER.py
			OSPM_RUN.py
		Observables
			OSPM_Observables_Stellar.py 
		Physics
			OSPM_Physics.py
			OSPM_PhysicsEngine.py
			OSPM_Physics_Spherical.jl
			OSPM_Physics_Sphericalv1.jl (archive it)
		Plotting
			plot_analysis.py
		Solvers
			OSPM_Solver_stellar.py
		_Archive_
		workspace.code-workspace 
	Utils
		Clear
		isitrunnin
		plot
		start
	_Setup
		Readme
		__init__.py
		requirements.txt

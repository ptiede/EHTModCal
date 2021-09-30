using EHTModelStacker


function makeh5val(dir)
  sdirs = filter(isdir, readdir(dir, join=true))
  println(sdirs)
  for sd in sdirs
    ssdirs = readdir(sd, join=true)
    for ssd in ssdirs
      println(ssd)
      try
        make_hdf5_chain_rose(ssd, joinpath(ssd, "chain.h5"))
      catch e
        println("Error on $ssd")
        println(e)
      end
    end
  end
end

function makeh5(dir)
  sdirs = filter(isdir, readdir(dir, join=true))
  println(sdirs)
  for sd in sdirs
    ssdirs = readdir(joinpath(sd, "noisefrac0.02"), join=true)
    for ssd in ssdirs
      println(ssd)
      try
        make_hdf5_chain_rose(ssd, joinpath(ssd, "chain.h5"))
      catch e
        println("Error on $ssd")
        println(e)
      end
    end
  end
end

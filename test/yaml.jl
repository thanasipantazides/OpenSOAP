using OpenSOAP
import YAML

function yaml_main(config::String)
    targets = load_mission(config)
    # println(targets)
end
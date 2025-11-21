
"""
Test script for the pyMBE DFManagement API
with hard-coded pmb_type in templates and instances.
"""

from pyMBE.storage.df_management import _DFManagement
from pyMBE.storage.templates.particle import ParticleTemplate
from pyMBE.storage.instances.particle import ParticleInstance


def main():
    print("=== Testing DFManagement API ===")

    pmb = _DFManagement()

    # -------------------------------------------------------
    # 1. Create two particle templates
    # -------------------------------------------------------
    tpl_A = ParticleTemplate(
        name="WCA_bead",
        sigma=1.0,
        epsilon=1.0
    )

    tpl_B = ParticleTemplate(
        name="LJ_bead",
        sigma=0.5,
        epsilon=2.0
    )

    # Register templates
    print("\nRegistering templates...")
    pmb.register_template(tpl_A)
    pmb.register_template(tpl_B)

    print("Registered templates:")
    for tname in pmb.templates["particle"]:
        print(f" - {tname}")    

    # -------------------------------------------------------
    # 2. Create particle instances
    # -------------------------------------------------------
    inst_1 = ParticleInstance(
        particle_id=1,  
        name="WCA_bead"
    )   
    inst_2 = ParticleInstance(
        particle_id=2,  
        name="LJ_bead"
    )

    inst_3 = ParticleInstance(
        particle_id=3,  
        name="LJ_bead"
    )

    # Register instances
    print("\nRegistering instances...")
    pmb.register_instance(inst_1)
    pmb.register_instance(inst_2)
    pmb.register_instance(inst_3)
    print("Registered instances:")
    print(pmb.instances)
    
if __name__ == "__main__":
    main()
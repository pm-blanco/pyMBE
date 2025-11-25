# test.py
from pyMBE.storage.df_management import _DFManagement
from pyMBE.storage.templates.particle import ParticleTemplate, ParticleState
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant

def main():
    db = _DFManagement()


    # ============================================================
    # 1. CREATE PARTICLE TEMPLATES + STATES
    # ============================================================

    # A particle (acid)
    tpl_A = ParticleTemplate(name="A", sigma=3.5, epsilon=0.2)
    tpl_A.add_state(ParticleState(name="HA", charge=0, es_type=0))
    tpl_A.add_state(ParticleState(name="A-", charge=-1, es_type=1))

    # H+ particle (single-state)
    tpl_H = ParticleTemplate(name="H", sigma=3.5, epsilon=0.2)
    tpl_H.add_state(ParticleState(name="H+", charge=+1, es_type=2))

    # Register templates
    db.register_template(tpl_A)
    db.register_template(tpl_H)

    # ============================================================
    # 2. CREATE INSTANCES (optional for testing)
    # ============================================================

    inst1 = ParticleInstance(name="A", particle_id=1, state_name="HA")
    inst2 = ParticleInstance(name="A", particle_id=2, state_name="A-")
    inst3 = ParticleInstance(name="H", particle_id=3, state_name="H+")

    db.register_instance(inst1)
    db.register_instance(inst2)
    db.register_instance(inst3)

    # ============================================================
    # 3. DEFINE A REACTION:  HA <-> A- + H+
    # ============================================================

    rx = Reaction(
        name="acid_dissociation",
        constant=1e-5,
        reaction_type="acid/base",
        participants=[
            ReactionParticipant(particle_name="A", state_name="HA", coefficient=-1),
            ReactionParticipant(particle_name="A", state_name="A-", coefficient=+1),
            ReactionParticipant(particle_name="H", state_name="H+", coefficient=+1),
        ],
    )

    db.register_reaction(rx)

    # ============================================================
    # 4. PRINT DATAFRAMES
    # ============================================================

    print("\n=== Templates DataFrame ===")
    print(db.get_templates_df())

    print("\n=== Instances DataFrame ===")
    print(db.get_instances_df())

    print("\n=== Reactions DataFrame ===")
    print(db.get_reactions_df())

if __name__ == "__main__":
    main()


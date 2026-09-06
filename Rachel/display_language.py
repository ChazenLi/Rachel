"""English presentation text for stable Rachel runtime identifiers."""

from __future__ import annotations

from typing import Any, Dict


DEFAULT_ENGLISH_SELF_PROMPT = (
    "Use English for all user-facing planning prose and explanations unless "
    "the user explicitly requests another language. Preserve SMILES, schema "
    "identifiers, enum values, and exact source quotations unchanged."
)


_EXPERIENCE_CARD_TEXT: Dict[str, Dict[str, str]] = {
    "exp_paal_knorr_deep_ring_warning": {
        "one_line": (
            "Deep Paal-Knorr disconnections of mature heteroaromatic scaffolds "
            "are high-risk by default and often signal a template-driven false "
            "deep disconnection."
        ),
        "action_prompt": (
            "First ask whether this heterocycle should be retained as an "
            "advanced terminal, then decide whether to test a Paal-Knorr route."
        ),
        "avoid": (
            "Do not treat `[NH2]` or an incomplete small-molecule representation "
            "as a high-quality terminal."
        ),
    },
    "exp_snar_electron_poor_heteroaryl": {
        "one_line": (
            "When an electron-poor fluoroheteroarene matches a suitable "
            "N-nucleophile, SNAr may be preferred over metal-catalyzed C-N coupling."
        ),
        "action_prompt": (
            "Check the leaving group, ortho/para activation, steric demand of the "
            "N-nucleophile, and competing nucleophilic sites."
        ),
        "avoid": (
            "Do not mechanically switch to Buchwald, Ullmann, or Chan-Lam coupling "
            "when SNAr is an evident match."
        ),
    },
    "exp_organometallic_precursor_normalization": {
        "one_line": (
            "Keep the organometallic reagent in the current reaction; its source "
            "hint should become a separate real upstream step rather than replace "
            "the current reagent."
        ),
        "action_prompt": (
            "When precursor_normalization appears, keep current_reagent in the "
            "current reaction and treat upstream_source_precursors as source "
            "candidates; compare direct metal insertion, transmetalation, or another "
            "realistic one-step preparation before deciding whether to create a "
            "continuation."
        ),
        "avoid": (
            "Do not replace the C-metal reagent required by the current reaction "
            "with an organic halide plus metal and then commit the same step; do not "
            "treat a source hint as the only mechanism, and do not treat a reactive "
            "organometallic as a final simple starting material."
        ),
    },
    "exp_protection_is_tree_node": {
        "one_line": (
            "Protection and deprotection are real route nodes, not reasoning annotations."
        ),
        "action_prompt": (
            "If conditions conflict with an unprotected functional group, represent "
            "protection or deprotection as an independent action and route step."
        ),
        "avoid": (
            "Do not write 'assume protected' in commit reasoning without adding the "
            "event to the tree."
        ),
    },
    "exp_route_plan_revision": {
        "one_line": (
            "When a better route emerges or evidence conflicts with the existing "
            "plan, make a concise route_plan revision instead of drifting silently."
        ),
        "action_prompt": (
            "Use route_plan(...) to record revision_reason and update the thesis, key "
            "disconnections, and structures to preserve before resuming the standard "
            "action/sandbox loop."
        ),
        "avoid": (
            "Do not mix the old plan and the new choice in reasoning without updating state."
        ),
    },
    "exp_chemist_guidance_is_direction_not_evidence": {
        "one_line": (
            "Chemist guidance is high-priority route direction, not validation evidence."
        ),
        "action_prompt": (
            "Convert natural-language guidance into site/reaction/precursor constraints; "
            "when no system action matches, use propose_action and then try_action."
        ),
        "avoid": (
            "Do not skip candidate comparison, sandboxing, or commit audit because an "
            "expert specified the direction."
        ),
    },
}


_DANGEROUS_COMBO_TEXT = {
    "nitro_amine_explosive": (
        "Nitro and amine combination: potential explosive hazard; tightly control "
        "reaction conditions."
    ),
    "azide_carbonyl_curtius": (
        "Azide and carbonyl combination: Curtius rearrangement hazard; apply "
        "appropriate safety controls."
    ),
    "epoxide_amine_opening": (
        "Epoxide and amine present together: risk of uncontrolled ring opening."
    ),
    "epoxide_thiol_opening": (
        "Epoxide and thiol present together: risk of rapid ring opening."
    ),
    "acyl_halide_moisture": (
        "Acyl halide and hydroxyl group present together: high reactivity and "
        "possible unintended acylation."
    ),
    "anhydride_nucleophile": (
        "Anhydride and amine present together: competing acylation risk."
    ),
    "diazo_heat_sensitive": (
        "Diazo and nitro groups present together: extreme instability and potential "
        "explosion hazard."
    ),
    "peroxide_reducing_agent": (
        "Peroxide and thiol present together: risk of a vigorous redox reaction."
    ),
    "organolithium_protic": (
        "Organolithium reagent and proton source present together: risk of a violent "
        "exothermic quench; strict anhydrous conditions are required."
    ),
    "diazo_strong_acid": (
        "Diazo compound with strong acid or nitro functionality: thermally unstable "
        "with potential explosion hazard (Barton/Shapiro intermediate)."
    ),
    "tempo_strong_reductant": (
        "TEMPO or another nitroxide radical with a strong reductant such as an "
        "organolithium reagent: radical quenching and reaction failure."
    ),
}


_SELECTIVITY_CONFLICT_TEXT = {
    ("aldehyde", "primary_amine", "imine_formation"): (
        "Aldehydes and primary amines spontaneously condense to imines (Schiff "
        "bases); protect one partner or tightly control pH and temperature."
    ),
    ("aldehyde", "secondary_amine", "enamine_formation"): (
        "Aldehydes and secondary amines may form enamines; control conditions or "
        "protect the aldehyde."
    ),
    ("aldehyde", "thiol", "hemithioacetal_formation"): (
        "Aldehydes and thiols rapidly form thermodynamically favored hemithioacetals; "
        "protect one partner."
    ),
    ("aldehyde", "ketone", "electrophilicity_competition"): (
        "Aldehydes are about 100 times more electrophilic than ketones; address "
        "selectivity in Grignard, reduction, and Wittig reactions."
    ),
    ("epoxide", "primary_amine", "ring_opening"): (
        "Epoxides and primary amines can undergo uncontrolled SN2 ring opening; "
        "protect one partner."
    ),
    ("epoxide", "secondary_amine", "ring_opening"): (
        "Epoxides and secondary amines can undergo SN2 ring opening; protect one partner."
    ),
    ("epoxide", "thiol", "ring_opening"): (
        "Epoxides and thiols undergo rapid ring opening because thiols are strongly "
        "nucleophilic; protect one partner."
    ),
    ("epoxide", "secondary_alcohol", "ring_opening"): (
        "Alcohols can open epoxides under acidic or basic conditions; review the "
        "reaction conditions."
    ),
    ("epoxide", "carboxylic_acid", "ring_opening"): (
        "Carboxylic acids can promote epoxide opening to hydroxy esters; control the "
        "reaction conditions."
    ),
    ("acyl_halide", "primary_amine", "acylation_competition"): (
        "Acyl halides react rapidly and irreversibly with amines; protect the amine "
        "when this is not the intended reaction."
    ),
    ("acyl_halide", "primary_alcohol", "acylation_competition"): (
        "Acyl halides rapidly esterify alcohols; protect the alcohol when this is not "
        "the intended reaction."
    ),
    ("carboxylic_acid", "primary_alcohol", "esterification"): (
        "Carboxylic acids and alcohols can esterify under activating conditions; "
        "control the activating reagent."
    ),
    ("carboxylic_acid", "secondary_alcohol", "esterification"): (
        "Carboxylic acids and secondary alcohols can esterify under activating "
        "conditions; control the activation."
    ),
    ("carboxylic_acid", "phenol", "esterification"): (
        "Carboxylic acids and phenols can esterify under activating conditions; "
        "control the activation."
    ),
    ("carboxylic_acid", "primary_amine", "salt_or_amidation"): (
        "Carboxylic acids and primary amines form ammonium salts and may undergo "
        "direct amidation after activation; control the activator or protect one partner."
    ),
    ("michael_acceptor", "primary_amine", "conjugate_addition"): (
        "Michael acceptors and primary amines can undergo 1,4-conjugate addition; "
        "protect the amine or control stoichiometry."
    ),
    ("michael_acceptor", "thiol", "thiol_michael"): (
        "Thiol-Michael addition can be extremely fast (about 10^3 M^-1 s^-1); "
        "protect one partner."
    ),
    ("michael_acceptor", "secondary_amine", "conjugate_addition"): (
        "Michael acceptors and secondary amines may undergo conjugate addition."
    ),
    ("alkene", "alkyne", "selectivity_competition"): (
        "When alkenes and alkynes coexist, hydrogenation and oxidation selectivity "
        "requires control; a Lindlar catalyst can selectively reduce an alkyne."
    ),
    ("alkene", "aldehyde", "conjugate_addition"): (
        "Conjugated enals may undergo unintended 1,4-addition."
    ),
    ("halide_aryl", "primary_amine", "buchwald_side_reaction"): (
        "A free amine may participate in side reactions under Buchwald coupling "
        "conditions; consider protection."
    ),
    ("halide_aryl", "boronic_acid", "protodeboronation_risk"): (
        "When an aryl halide and boronic acid coexist, the boronic acid may undergo "
        "deboronation; control the reaction order."
    ),
    ("primary_alcohol", "primary_amine", "nucleophile_competition"): (
        "Alcohols and amines are both nucleophiles; amine acylation is usually about "
        "10^3 times faster, but selectivity must be confirmed."
    ),
    ("phenol", "primary_amine", "acylation_competition"): (
        "Phenols and amines compete in acylation; amine acylation is kinetically "
        "favored, but thermodynamic O-acylation of the phenol is not negligible."
    ),
    ("halide_alkyl", "primary_alcohol", "sn2_competition"): (
        "SN2 reaction between an alkyl halide and alcohol may compete with the "
        "intended transformation."
    ),
    ("nitrile", "primary_amine", "amidine_formation"): (
        "Nitriles and amines may form amidines under strong acid or Lewis acid conditions."
    ),
    ("alkyne", "primary_amine", "hydroamination"): (
        "Alkynes and amines may undergo hydroamination under transition-metal catalysis."
    ),
    ("ester", "primary_amine", "aminolysis"): (
        "Esters and amines can undergo aminolysis to form amides, especially with "
        "heating or under basic conditions."
    ),
    ("thiol", "disulfide", "thiol_disulfide_exchange"): (
        "Free thiols undergo exchange with disulfides; protect the thiol or control "
        "the redox environment."
    ),
    ("thiol", "halide_alkyl", "alkylation"): (
        "Strongly nucleophilic thiols may undergo unintended alkylation with alkyl halides."
    ),
    ("thiol", "primary_alcohol", "nucleophile_competition"): (
        "Thiols are about 10^3 times more nucleophilic than alcohols and can dominate "
        "competing nucleophilic reactions."
    ),
    ("azide", "alkyne", "click_reaction"): (
        "Azides and alkynes undergo Cu(I)-catalyzed Huisgen cycloaddition; account for "
        "this when it is not the intended reaction."
    ),
    ("sulfonate_ester", "primary_amine", "nucleophilic_displacement"): (
        "Sulfonate esters are good leaving groups and amines can displace them directly; "
        "protect the amine when needed."
    ),
    ("sulfonate_ester", "thiol", "nucleophilic_displacement"): (
        "Sulfonate esters undergo rapid nucleophilic displacement by thiols."
    ),
}


def project_experience_card(card: Dict[str, Any]) -> Dict[str, Any]:
    projected = dict(card)
    projected.update(_EXPERIENCE_CARD_TEXT.get(str(card.get("id", "")), {}))
    return projected


def dangerous_combo_text(combo_key: str, fallback: str = "") -> str:
    return _DANGEROUS_COMBO_TEXT.get(str(combo_key), str(fallback or ""))


def selectivity_conflict_text(
    group_a: str,
    group_b: str,
    reaction_type: str,
    fallback: str = "",
) -> str:
    key = (str(group_a), str(group_b), str(reaction_type))
    return _SELECTIVITY_CONFLICT_TEXT.get(key, str(fallback or ""))

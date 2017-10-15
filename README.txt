HROCH: Heuristic reconstruction of cluster histories.

Pôvodná verzia:
- https://github.com/jablkoj/hroch
- https://people.ksp.sk/~janoh/diplomovka.php

Skompilujeme príkazom make.

Odporúčané použitie:
./hroch.bin --solve SÚBOR_ATÓMOV PRIEČINOK_STROMOV POČET [STRATÉGIA]

Zrekonštruuje daný POČET histórií sekvencie danej súborom SÚBOR_ATÓMOV.
Vypíše štatistiky a množinu najkratších histórií.
V priešinku outputs vytvorí nový súbor kde vypíše všetky histórie.
Podporované sú nasledovné stratégie:
    1: bez čerešňovitosti a skórovania
    2: s čerešňovitosťou ale bez skórovania
    5: bez čerešňovitosti, jednoduché skórovanie
    6: s čerešňovitosťou, jednoduché skórovanie
    7: pokročílé skórovanie s čerešňovitosťou
    8: [predvolené] pokročílé striktné skórovanie s čerešňovitosťou


Program sa dá použiť aj na mnoho ďalších nezdokumentovaných vecí,
odporúčame konzultovať s autorom.

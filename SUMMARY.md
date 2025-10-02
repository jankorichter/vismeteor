# vismeteor – Kurzüberblick

## Zweck und Funktionen
`vismeteor` bündelt Analysewerkzeuge für visuelle Meteorbeobachtungen aus der Visual Meteor Database (VMDB) der International Meteor Organization. Neben Standardauswertungen unterstützt das Paket u. a. die Parameterbestimmung für geometrische und ideale Modelle sowie die Aufbereitung eigener Beobachtungsdaten.

## Datenquellen
Standardmäßig greift das Paket auf Datensätze aus dem begleitenden Projekt `imo-vmdb` zurück. Die Funktionen sind jedoch so gestaltet, dass auch alternative oder eigene Datenquellen eingebunden werden können.

## Installation und Einstieg
Die stabile Version ist auf CRAN verfügbar und lässt sich mit `install.packages("vismeteor")` installieren. Für den Einstieg stehen mehrere Vignetten bereit, darunter `vignette("vismeteor")` als Gesamtüberblick sowie `vignette("vmgeom")` und `vignette("vmideal")` zur Parameterschätzung spezieller Modelle.

## Mitwirkung, Lizenz und Kontakt
Beiträge der Community sind willkommen. Der empfohlene Workflow umfasst Fork, Feature-Branch, Commit, Push und Pull Request. Das Paket steht unter der MIT-Lizenz; Ansprechpartner ist Janko Richter (janko@richtej.de).

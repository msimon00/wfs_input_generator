#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DESCRIPTION

:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2013
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""
from obspy import readEvents
from obspy.core import AttribDict
from obspy.core.event import Event
from obspy.xseed import Parser


class Input_File_Generator(object):
    def __init__(self):
        self.config = AttribDict()
        self._events = []
        self._stations = []

    def add_events(self, events):
        """
        Add one or more events to the input file generator. Most inversions
        should specify only one event but some codes can deal with multiple
        events.

        Can currently deal with QuakeML files and obspy.core.event.Event
        objects.

        :type events: list or obspy.core.event.Catalog object
        :param events: A list of filenames, a list of obspy.core.event.Event
            objects, or an obspy.core.event.Catalog object.
        """
        # Thin wrapper to be able to also treat single events or filenames.
        if isinstance(events, Event) or not hasattr(events, "__iter__"):
            events = list(events)

        for event in events:
            if isinstance(event, Event):
                self._parse_event(event)
                continue
            try:
                cat = readEvents(event)
            except:
                msg = "Could not read %s." % str(event)
                raise TypeError(msg)
            for event in cat:
                self._parse_event(event)

    def add_stations(self, stations):
        """
        Add the desired output stations to the input file generator.

        Can currently deal with SEED/XML-SEED files and dictionaries of the
        following form:

            {"latitude": 123.4,
             "longitude": 123.4,
             "elevation_in_m": 123.4,
             "local_depth_in_m": 123.4,
             "id": "network_code.station_code"}

        :type stations: List of filenames, list of dictionaries or a single
            filename, single dictionary.
        :param stations: The stations for which output files should be
            generated.
        """
        # Thin wrapper to enable single element treatment.
        if isinstance(stations, dict) or not hasattr(stations, "__iter__"):
            stations = list(stations)
        pass


    def write(self, events):
        pass

    def _parse_event(self, event):
        """
        Check and parse events.

        Each event at least needs to have an origin and a moment tensor,
        otherwise an error will be raised.
        """
        # Do a lot of checks first to be able to give descriptive error
        # messages.
        if not event.origins:
            msg = "Each event needs to have an origin."
            raise ValueError(msg)
        if not event.focal_mechanisms:
            msg = "Each event needs to have a focal mechanism."
            raise ValueError(msg)
        # Choose either the preferred origin or the first one.
        origin = event.preferred_origin() or event.origins[0]
        # Same with the focal mechanism.
        foc_mec = event.preferred_focal_mechanism() or \
            event.focal_mechanisms[0]
        # Origin needs to have latitude, longitude, depth and time
        if not origin.latitude or not origin.longitude or not origin.depth \
            or not origin.time:
            msg = ("Every event origin needs to have latitude, longitude, "
                "depth and time")
            raise ValueError(msg)
        # The focal mechanism of course needs to have a moment tensor.
        if not foc_mec.moment_tensor or not foc_mec.moment_tensor.tensor:
            msg = "Every event needs to have a moment tensor."
            raise ValueError(msg)
        # Also all six components need to be specified.
        mt = foc_mec.moment_tensor.tensor
        if not mt.m_rr or not mt.m_tt or not mt.m_pp or not mt.m_rt \
            or not mt.m_rp or not mt.m_tp:
            msg = "Every event needs all six moment tensor components."
            raise ValueError(msg)

        # Now the event should be valid.
        self.events.append({
            "latitude": origin.latitude,
            "longitude": origin.longitude,
            "depth_in_km": origin.depth,
            "origin_time": origin.time,
            "m_rr": mt.m_rr,
            "m_tt": mt.m_tt,
            "m_pp": mt.m_pp,
            "m_rt": mt.m_rt,
            "m_rp": mt.m_rp,
            "m_tp": mt.m_tp})